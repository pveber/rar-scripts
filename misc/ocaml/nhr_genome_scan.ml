open Batteries
open Biocaml
open Printf

let selected_motifs = List.concat [
  List.init 10 (fun i -> `direct, i) ;
  List.init 10 (fun i -> `everted, i) ;
  List.init 10 (fun i -> `inverted, i) ;
]

let selected_motifs = [ `direct, 5 ]

let scan_motif m seq =
  Nhre.matches_of_sequence [ m ] Nhre.estimated_counts 20. seq

let adapter item = Oregon.Fasta.(item.comment, item.seq)

let scan (chr, seq) =
  chr, 
  List.map
    (fun m -> scan_motif m seq)
    selected_motifs

let scana (chr, seq) =
  chr, 
  String.length seq

let _ = 
  (Fasta.enum_of_file B.Genome.sequence#path |> snd |> Enum.take 1) /@ (scan) 
  |> List.of_enum
    
  

