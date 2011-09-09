open Batteries
open Printf
open Oregon
open Sle


let ranking radius =
  let fa = Ucsc.fasta_of_bed `mm9 (Manual_chipseq_clustering.gerard_selected) in
  let ctrl = Fasta.shuffle fa in
  sprintf "drerir-%d.pdf" radius,
  Motif_library.rank Selected_motifs.drerir fa ctrl


let archive = 
  let r = 
    List.map ranking [ 50 ; 150 ; 250 ; 350 ; 450 ] 
    |> List.map (Tuple2.mapn identity (fun x -> x#abstract))
  in
  Archive.make r

