(** API to fetch/use data from ensembl *)

type species = [
| `mus_musculus
]

val fetch_gtf : release:int -> species:species -> path:string -> unit




















