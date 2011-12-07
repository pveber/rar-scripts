open Oregon
open Genome
open Oregon.Target

type sample = B.Rnaseq.sample
type condition = string * sample list * sample list

type annot = {
  id : string ;
  baseMeanA : float array ; (* for all conditions *)
  baseMeanB : float array ;
  padj      : float array ;
}

val make : 
  condition list -> Gene.t array value -> (condition array * annot array) value

val tsv : string -> (condition array * annot array) value -> unit
