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


let pval_of_string s = 
  try float_of_string s
  with _ -> 1.

let log2_of_string s = 
  try float_of_string s
  with _ -> 0.

let base_mean_of_string = log2_of_string

let parse_row r = {
  ensembl_gene_id = r.(0) ;
  gene_symbol = r.(1) ;
  chr = "chr" ^ r.(6) ;
  chr_start = int_of_string r.(7) ;
  chr_end = int_of_string r.(8) ;
  chr_strand = (match r.(9) with 
  | "1" -> `Sense
  | "-1" -> `Antisense
  | s -> print_endline s ; assert false) ;
  base_mean_6H = base_mean_of_string r.(11) ;
  base_mean_12H = base_mean_of_string r.(17) ;
  base_mean_24H = base_mean_of_string r.(22) ;
  base_mean_36H = base_mean_of_string r.(29) ;
  base_mean_48H = base_mean_of_string r.(32) ;
  log2_fold_6H = log2_of_string r.(13) ;
  log2_fold_12H = log2_of_string r.(19) ;
  log2_fold_24H = log2_of_string r.(24) ;
  log2_fold_36H = log2_of_string r.(29) ;
  log2_fold_48H = log2_of_string r.(34) ;
  pval_6H = pval_of_string r.(14) ;
  pval_12H = pval_of_string r.(20) ;
  pval_24H = pval_of_string r.(25) ;
  pval_36H = pval_of_string r.(30) ;
  pval_48H = pval_of_string r.(35) ;
  padj_6H = pval_of_string r.(15) ;
  padj_12H = pval_of_string r.(21) ;
  padj_24H = pval_of_string r.(26) ;
  padj_36H = pval_of_string r.(31) ;
  padj_48H = pval_of_string r.(36) ;
}

let loc_of_row r = Biocaml_range.(r.chr, make r.chr_start r.chr_end)

let of_file fn = 
  (File.lines_of fn |> Enum.skip 1)
  /@ Tsv.row_of_string
  /@ parse_row

let detection_level = 9.99e-9

let expressed g =
  (abs_float g.base_mean_6H >= detection_level) ||
  (abs_float g.base_mean_12H >= detection_level) ||
  (abs_float g.base_mean_24H >= detection_level) ||
  (abs_float g.base_mean_36H >= detection_level) ||
  (abs_float g.base_mean_48H >= detection_level)

let expressed_6H g =
  (abs_float g.base_mean_6H >= detection_level)

let expressed_12H g =
  (abs_float g.base_mean_6H >= detection_level) ||
  (abs_float g.base_mean_12H >= detection_level)

let expressed_24H g =
  (abs_float g.base_mean_6H >= detection_level) ||
  (abs_float g.base_mean_12H >= detection_level) ||
  (abs_float g.base_mean_24H >= detection_level)

let expressed_36H g =
  (abs_float g.base_mean_6H >= detection_level) ||
  (abs_float g.base_mean_12H >= detection_level) ||
  (abs_float g.base_mean_24H >= detection_level) ||
  (abs_float g.base_mean_36H >= detection_level)

let modulated g =
  (abs_float g.log2_fold_6H >=1. && g.padj_6H <= 0.05) ||
  (abs_float g.log2_fold_12H >=1. && g.padj_12H <= 0.05) ||
  (abs_float g.log2_fold_24H >=1. && g.padj_24H <= 0.05) ||
  (abs_float g.log2_fold_36H >=1. && g.padj_36H <= 0.05) ||
  (abs_float g.log2_fold_48H >=1. && g.padj_48H <= 0.05)

let upregulated g =
  (g.log2_fold_6H >=1. && g.padj_6H <= 0.05) ||
  (g.log2_fold_12H >=1. && g.padj_12H <= 0.05) ||
  (g.log2_fold_24H >=1. && g.padj_24H <= 0.05) ||
  (g.log2_fold_36H >=1. && g.padj_36H <= 0.05) ||
  (g.log2_fold_48H >=1. && g.padj_48H <= 0.05)

let upregulated_6H g =
  (g.log2_fold_6H >=1. && g.padj_6H <= 0.05)

let upregulated_12H g =
  (g.log2_fold_12H >=1. && g.padj_12H <= 0.05)

let upregulated_24H g =
  (g.log2_fold_24H >=1. && g.padj_24H <= 0.05)

let upregulated_36H g =
  (g.log2_fold_36H >=1. && g.padj_36H <= 0.05)

let upregulated_48H g =
  (g.log2_fold_48H >=1. && g.padj_48H <= 0.05)

let downregulated g =
  (g.log2_fold_6H <= -1. && g.padj_6H <= 0.05) ||
  (g.log2_fold_12H <= -1. && g.padj_12H <= 0.05) ||
  (g.log2_fold_24H <= -1. && g.padj_24H <= 0.05) ||
  (g.log2_fold_36H <= -1. && g.padj_36H <= 0.05) ||
  (g.log2_fold_48H <= -1. && g.padj_48H <= 0.05)

let downregulated_6H g =
  (g.log2_fold_6H <= -1. && g.padj_6H <= 0.05)

let downregulated_12H g =
  (g.log2_fold_12H <= -1. && g.padj_12H <= 0.05)

let downregulated_24H g =
  (g.log2_fold_24H <= -1. && g.padj_24H <= 0.05)

let downregulated_36H g =
  (g.log2_fold_36H <= -1. && g.padj_36H <= 0.05)

let downregulated_48H g =
  (g.log2_fold_48H <= -1. && g.padj_48H <= 0.05)




















