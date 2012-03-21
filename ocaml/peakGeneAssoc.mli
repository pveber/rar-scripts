(** 
    Association between chipseq regions and genes, distance histograms
*)
open Batteries
open Genome

val closest_ensembl_gene : Location.t Enum.t -> (Location.t * string * int) Enum.t

val panRAR_regions_to_gene_tsv : string -> unit




















