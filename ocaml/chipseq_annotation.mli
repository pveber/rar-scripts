(** annotation d'un ensemble de régions par les résultats de chipseq *)

open Oregon
open Genome
open Oregon.Target

type 'a annot = {
  loc : Location.t ;
  counts : ('a, int) HFun.t ;
  coverage : ('a, float) HFun.t ;
  pvalue : ('a, float) HFun.t
}

val make : 
  (B.Chipseq.sample * B.Chipseq.sample) list -> Location.t array value -> B.Chipseq.sample annot array value

val tsv : string -> B.Chipseq.sample annot array -> unit
