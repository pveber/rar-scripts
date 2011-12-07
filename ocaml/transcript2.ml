open Genome

type t = {
  id : string ;
  gene_symbol : string ;
  description : string ;
  exons : Location.t list ;
  strand : [ `plus | `minus ] ;
}
