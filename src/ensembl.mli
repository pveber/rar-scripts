open Batteries
open Biocaml
open GenomeMap
open Guizmin_bioinfo.MBSchema

(** API to fetch/use data from ensembl *)

type species = [
| `mus_musculus
]

val fetch_gtf : release:int -> species:species -> path:string -> unit

val transcripts_enum : Gtf.item Enum.t -> Transcript.t Enum.t
(** annotation is gene_id, transcript_id *)

val tss_map_of_gtf : Gtf.item Enum.t -> (string, Transcript.t) LMap.t
(** annotation is gene_id, transcript_id *)

val genes_enum : Gtf.item Enum.t -> Gene.t Enum.t




















