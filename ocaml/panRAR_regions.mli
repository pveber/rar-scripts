(** PanRAR selected regions as defined with the island algorithm,
    applied on panRAR ChIP-seq data *)
open Genome
open Oregon.Target

val conditions : B.Chipseq.sample list
val all_of_them : Location.t array value
val gerard_selection : Location.t array value




















