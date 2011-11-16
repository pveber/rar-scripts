(* Test qui ne marche pas du tout, candidat à la suppression ... *)
open Batteries
open Printf
open Genome
open Oregon
open Sle

let ball n loc = Location.(
  let i = (loc.st + loc.ed) / 2 
  and w = (loc.ed - loc.st) / 2 in
  make loc.chr (i - w) (i + w)
)

let ranking radius =
  let path = Filename.temp_file "guiz" (string_of_int (Random.int 1000000)) in
  let locz = 
    Tsv.enum Manual_chipseq_clustering.gerard_selected 
    /@ (fun x -> x#loc)
    /@ (ball radius) 
  in
  locz /@ Bed.unparser /@ Tsv.string_of_line |> File.write_lines path ;	(* UGLY ! *)		
  let bed = Target.F.input path Tsv.({ has_header = false ; parse = new Bed.base_parser }) in 
  let fa = Ucsc.fasta_of_bed `mm9 bed in
  let ctrl = Fasta.shuffle fa in
  sprintf "drerir-%d.pdf" radius,
  Motif_library.rank Selected_motifs.drerir fa ctrl


let archive = 
  let r = 
    List.map ranking [ 50 ; 150 ; 250 ; 350 ; 450 ] 
    |> List.map (Tuple2.mapn identity (fun x -> x#abstract))
  in
  Archive.make r

