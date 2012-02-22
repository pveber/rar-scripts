open Batteries
open Printf
open Biocaml
open Misc

let selected_motifs = 
  [ `direct, 0 ; `direct, 1 ; `direct, 2 ; `direct, 5 ] ::
  ((List.init 10 (fun i -> [ `direct, i ])) @
   (List.init 10 (fun i -> [`everted, i]))    @
   (List.init 10 (fun i -> [`inverted, i])))


let fa_input = Sys.argv.(1)
let output = Sys.argv.(2)
let _ = ignore (Sys.command ("mkdir -p " ^ output))

let fa =
  Fasta.enum_of_file fa_input
  |> snd
  |> List.of_enum

let sequences = List.map snd fa
let control_sequences = List.map string_shuffle sequences
let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) (*string_shuffle*) sequences

let roc_curve hexamer_counts selected_motifs =
  let matches = 
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) sequences
  and control_matches = 
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) control_sequences in
  let extract_best l = 
    List.fold_left (fun accu (_,_,_,x) -> max accu x) neg_infinity l
  in
  let pos = List.map extract_best matches
  and neg = List.map extract_best control_matches in
  Roc.(
    make ~pos:(List.enum pos) ~neg:(List.enum neg)
    /@ (fun (_,cm) -> false_positive_rate cm, sensitivity cm)
    |> List.of_enum
  ) 
  

let string_of_motifs ms =
  List.map Nhre.string_of_motif ms |> String.concat "-" 

let roc_curve_figure dir hexamer_counts selected_motifs = 
  let id = string_of_motifs selected_motifs in
  let path = dir ^ "/" ^ id ^ ".pdf" in
  let curve = roc_curve hexamer_counts selected_motifs in
  let auc = Roc.auc (List.enum curve) in
  let x = Array.of_list (List.map fst curve) 
  and y = Array.of_list (List.map snd curve) 
  and rp = R.make () in
  R.v rp "x" x ;
  R.v rp "y" y ;
  R.c rp "pdf('%s')" path ;
  R.c rp "plot(x,y,main = '%s occurrences, treatment VS control', xlab = 'Frequency in control sequences', ylab = 'Frequency in treatment sequences', type = 'l')" id ;
  R.c rp "abline(0,1,lty = 2)" ;
  R.c rp "dev.off()" ;
  R.close rp ;
  id, auc


let occurrences_tsv hexamer_counts = 
  let foreach selected_motifs = 
    let matches = 
      List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 4.) sequences
    in
    List.mapi 
      (fun i hits -> 
	List.enum hits
	/@ (fun (motif,sense,pos,score) -> sprintf "%d\t%s\t%s\t%d\t%.3f" i (Nhre.string_of_motif motif) (if sense = `sense then "+" else "-") pos score))
      matches
    |> List.enum
    |> Enum.concat
  in
  (List.enum (List.tl selected_motifs) /@ foreach)
  |> Enum.concat
  |> File.write_lines (output ^ "/predictions.tsv")


let res = 
  List.map
    (roc_curve_figure output Nhre.estimated_counts)
    selected_motifs

let () = 
  List.enum res
  /@ (fun (id, auc) -> sprintf "%s\t%f" id auc)
  |> File.write_lines (output ^ "/auc.tsv")

let () = occurrences_tsv Nhre.estimated_counts
