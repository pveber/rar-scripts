(** amélioration du module Transcript défini dans Genome qui est insuffisant *)
open Genome

type t = {
  id : string ;
  gene_symbol : string ;
  description : string ;
  exons : Location.t list ;
  strand : [ `plus | `minus ] ;
}
