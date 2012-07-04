open Guizmin

type output

val app : output dir pipeline

val rel_pos_bregions_fig : output dir pipeline -> [`pdf] file pipeline
(** relative position of binding regions with respect to the closest
    TSS *)
