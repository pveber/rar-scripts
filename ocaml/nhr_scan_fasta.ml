open Batteries
open Biocaml
open Printf

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

let print_hits hits = 
  let fields = 
    List.map
      (fun (motif, sense, i, s) -> 
	sprintf "%f" s)
      hits
  in
  String.concat "\t" fields
  |> print_endline

let () = 
  List.map Nhre.string_of_motif selected_motifs
  |> String.concat "\t"
  |> print_endline 

let () = 
  (snd (Fasta.enum_of_file Sys.argv.(1)) /@ scan)
  |> Enum.iter print_hits

