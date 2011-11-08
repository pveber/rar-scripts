open Oregon
open Genome
open Oregon.Target

type 'a annot = {
  coverage : ('a, int) HFun.t ;
  rpkm : ('a, float) HFun.t ;
  pvalue : ('a, float) HFun.t
}

val make : 
  (B.Chipseq.sample * B.Chipseq.sample) list -> Location.t array value -> B.Chipseq.sample annot array value

