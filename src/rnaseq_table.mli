(** 
    Parsing du tableau de synthèse d'Amandine pour le RNA-seq
*)

open Batteries

type row = {
  ensembl_gene_id : string ;
  gene_symbol : string ;
  chr : string ;
  chr_start : int ;
  chr_end : int ;
  chr_strand : [`Sense | `Antisense] ;
  base_mean_6H : float ;
  base_mean_12H : float ;
  base_mean_24H : float ;
  base_mean_36H : float ;
  base_mean_48H : float ;
  log2_fold_6H : float ;
  log2_fold_12H : float ;
  log2_fold_24H : float ;
  log2_fold_36H : float ;
  log2_fold_48H : float ;
  pval_6H : float ;
  pval_12H : float ;
  pval_24H : float ;
  pval_36H : float ;
  pval_48H : float ;
  padj_6H : float ;
  padj_12H : float ;
  padj_24H : float ;
  padj_36H : float ;
  padj_48H : float ;
}

val of_file : string -> row Enum.t

val loc_of_row : row -> string * Biocaml_range.t
(*val tss_loc_of_row : row -> string * Biocaml_range.t*)

val modulated : row -> bool
val expressed : row -> bool
val expressed_6H : row -> bool
val expressed_12H : row -> bool
val expressed_24H : row -> bool
val expressed_36H : row -> bool
val upregulated : row -> bool
val upregulated_6H : row -> bool
val upregulated_12H : row -> bool
val upregulated_24H : row -> bool
val upregulated_36H : row -> bool
val upregulated_48H : row -> bool
val downregulated : row -> bool
val downregulated_6H : row -> bool
val downregulated_12H : row -> bool
val downregulated_24H : row -> bool
val downregulated_36H : row -> bool
val downregulated_48H : row -> bool




















