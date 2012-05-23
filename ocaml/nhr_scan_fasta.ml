open Batteries
open Biocaml
open Printf

let fasta_file = Sys.argv.(1)
let output_dir = Sys.argv.(2)
let _ = ignore (Sys.command ("mkdir -p " ^ output_dir))

let selected_motifs = List.concat [
  List.init 10 (fun i -> `direct, i) ;
  List.init 10 (fun i -> `everted, i) ;
  List.init 10 (fun i -> `inverted, i) ;
]


let scan_motif m seq =
  let matches = Nhre.matches_of_sequence [ m ] Nhre.estimated_counts neg_infinity seq in
  let max ((_,_,_,sx) as x) ((_,_,_,sy) as y) = if sx > sy then x else y in
  List.reduce max matches

let scan (_,seq) =
  List.map
    (fun m -> scan_motif m seq)
    selected_motifs

let serialize_hits hits = 
  let fields = 
    List.map
      (fun (motif, sense, i, s) -> 
	sprintf "%f" s)
      hits
  in
  String.concat "\t" fields

let () = 
  (snd (Fasta.enum_of_file fasta_file) /@ scan)
  /@ serialize_hits
  |> Enum.append (Enum.singleton (List.map Nhre.string_of_motif selected_motifs |> String.concat "\t"))
  |> File.write_lines (output_dir ^ "/best_score.tsv")


let occurrences_tsv hexamer_counts = 
  let sequences = (snd (Fasta.enum_of_file fasta_file)) /@ snd |> List.of_enum in
  let foreach motif = 
    let matches = 
      List.map (Nhre.matches_of_sequence [motif] hexamer_counts 4.) sequences
    in
    List.mapi 
      (fun i hits -> 
	List.enum hits
	/@ (fun (motif,sense,pos,score) -> sprintf "%d\t%s\t%s\t%d\t%.3f" i (Nhre.string_of_motif motif) (if sense = `sense then "+" else "-") pos score))
      matches
    |> List.enum
    |> Enum.concat
  in
  (List.enum selected_motifs /@ foreach)
  |> Enum.concat
  |> File.write_lines (output_dir ^ "/predictions.tsv")

let () = occurrences_tsv Nhre.estimated_counts
















