open Batteries
open Genome
open Oregon
open Printf
open Scanf


let all_regions_fa = 
  let path = "manual/chipseq/regions/all_regions.bed" in
  (all_regions ()) /@ Bed.unparser /@ Tsv.string_of_line |> File.write_lines path ;	(* UGLY ! *)		
  let bed = Target.F.input path Tsv.({ has_header = false ; parse = new Bed.base_parser }) in 
  Ucsc.fasta_of_bed `mm9 bed

let make locz prefix = 
