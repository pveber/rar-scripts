open Guizmin

type output

val app : output dir pipeline

val fig_rel_pos_bregions : output dir pipeline -> [`pdf] file pipeline
(** relative position of binding regions with respect to the closest
    TSS *)

val fig_rel_pos_closest_bregion : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding region *)

val fig_expr_class_dist_wrt_closest_bregion_pos : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding regions *)

val fig_expr_class_dist_wrt_closest_bregion_pos_6H : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding regions *)

val fig_expr_class_dist_wrt_closest_bregion_pos_12H : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding regions *)

val fig_expr_class_dist_wrt_closest_bregion_pos_24H : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding regions *)

val fig_expr_class_dist_wrt_closest_bregion_pos_36H : output dir pipeline -> [`pdf] file pipeline
(** for all genes, relative position of the closest binding regions *)

val fig_upregulated_prop_wrt_closest_bregion_pos : output dir pipeline -> [`pdf] file pipeline

val fig_downregulated_prop_wrt_closest_bregion_pos : output dir pipeline -> [`pdf] file pipeline





















