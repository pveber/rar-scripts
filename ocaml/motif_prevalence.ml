open Batteries
open Printf
open Biocaml
open Misc

let selected_motifs = 
  [ `direct, 0 ; `direct, 1 ; `direct, 2 ; `direct, 5 ] ::
  ((List.init 10 (fun i -> [ `direct, i ])) @
   (List.init 10 (fun i -> [`everted, i]))    @
   (List.init 10 (fun i -> [`inverted, i])))


let output = "results/chipseq/regions/PanRAR_regions_motif_prevalence"
let _ = ignore (Sys.command ("mkdir -p " ^ output))

let rec mask m l = match m, l with 
| [], [] -> []
| true :: t, h' :: t' -> h' :: (mask t t')
| false :: t, h' :: t' -> mask t t'
| _ -> assert false

module CA = Chipseq_annotation

let chipseq_design = 
  List.concat [
    List.map (fun x -> x, `Input_1) PanRAR_regions.conditions ;
    [`F9_ATRA_panRAR_1,`F9_WT_panRAR_1; `F9_ATRA24_panRAR_1,`F9_WT_panRAR_1;`F9_ATRA48_panRAR_1,`F9_WT_panRAR_1] ;
    [`F9_ATRA2_panRXR_1,`F9_WT_panRXR_1;`F9_ATRA24_panRXR_1,`F9_WT_panRXR_1;`F9_ATRA48_panRXR_1,`F9_WT_panRXR_1] ;
  ]

let chipseq_annotation = (CA.make chipseq_design PanRAR_regions.gerard_selection)#value
let regions = Tuple3.first chipseq_annotation |> Array.to_list
let i2 = CA.i2 chipseq_annotation

let sequences = 
  PanRAR_regions.gerard_selection#value
  |> Oregon.Ucsc.genomic_sequences `mm9
  |> Array.to_list

let control_sequences = List.map string_shuffle sequences
let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) (*string_shuffle*) sequences

let roc_curve hexamer_counts selected_regions selected_motifs =
  let matches = 
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) (mask selected_regions sequences)
  and control_matches = 
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) (mask selected_regions control_sequences) in
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

let roc_curve_figure dir hexamer_counts selected_region_label selected_regions selected_motifs = 
  let id = string_of_motifs selected_motifs in
  let path = dir ^ "/" ^ selected_region_label ^ "." ^ id ^ ".pdf" in
  let curve = roc_curve hexamer_counts selected_regions selected_motifs in
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


let res_for_subset subset_label subset = 
  List.map
    (roc_curve_figure output Nhre.estimated_counts subset_label subset)
    selected_motifs
  |> List.enum
  |> Enum.map (fun (id, auc) -> sprintf "%s\t%f" id auc)
  |> File.write_lines (sprintf "%s/%s.auc.tsv" output subset_label)

let () = res_for_subset "all" (List.map (fun _ -> true) regions)
let () = 
  let filter ann = 
    ann.CA.pvalue.(i2 `F9_WT_panRAR_1 `Input_1) >= 7. &&
    ann.CA.pvalue.(i2 `F9_WT_panRXR_1 `Input_1) >= 7.
  in
  res_for_subset "T0" (List.map filter regions)

let () = 
  let filter ann = 
    ann.CA.pvalue.(i2 `F9_ATRA_panRAR_1 `Input_1) >= 7. &&
    ann.CA.pvalue.(i2 `F9_ATRA2_panRXR_1 `Input_1) >= 7.
  in
  res_for_subset "T2" (List.map filter regions)

let () = 
  let filter ann = 
    ann.CA.pvalue.(i2 `F9_ATRA24_panRAR_1 `Input_1) >= 7. &&
    ann.CA.pvalue.(i2 `F9_ATRA24_panRXR_1 `Input_1) >= 7.
  in
  res_for_subset "T24" (List.map filter regions)

let () = 
  let filter ann = 
    ann.CA.pvalue.(i2 `F9_ATRA48_panRAR_1 `Input_1) >= 7. &&
    ann.CA.pvalue.(i2 `F9_ATRA48_panRXR_1 `Input_1) >= 7.
  in
  res_for_subset "T48" (List.map filter regions)

let () = occurrences_tsv Nhre.estimated_counts


















