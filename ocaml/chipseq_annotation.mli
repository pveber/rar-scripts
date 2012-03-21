(** annotation d'un ensemble de régions par les résultats de chipseq *)

open Oregon
open Genome
open Oregon.Target

type annot = {
  loc : Location.t ;
  counts : int array ;
  coverage : float array ; 
  pvalue : float array ;
}

type t = annot array * B.Chipseq.sample array * (B.Chipseq.sample * B.Chipseq.sample) array

val make : 
  (B.Chipseq.sample * B.Chipseq.sample) list -> Location.t array value -> t value

val tsv : string -> t -> unit

val i :  t -> B.Chipseq.sample -> int
val i2 : t -> B.Chipseq.sample -> B.Chipseq.sample -> int


















