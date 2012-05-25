open Batteries
open Biocaml
open GenomeMap

(** API to fetch/use data from ensembl *)

type species = [
| `mus_musculus
]

val fetch_gtf : release:int -> species:species -> path:string -> unit

val tss_map_of_gtf : Gtf.item Enum.t -> (string, string) LMap.t




















