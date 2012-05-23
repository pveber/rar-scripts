(** 
    Association between chipseq regions and genes, distance histograms
*)
open Batteries
open Genome

val closest_ensembl_gene : Location.t Enum.t -> (Location.t * string * int * int) Enum.t
(** maps locations to the gene that has the closest TSS, giving its
    location, gene_id, distance in bp and its rank *)

val panRAR_regions_to_gene_tsv : string -> unit




















