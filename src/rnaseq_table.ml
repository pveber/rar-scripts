open Batteries

type row = {
  ensembl_gene_id : string ;
  gene_symbol : string ;
  chr : string ;
  chr_start : int ;
  chr_end : int ;
  chr_strand : [`Sense | `Antisense] ;
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

let tss_loc_of_row r = Biocaml_range.(
  let pos = match r.chr_strand with
  | `Sense -> r.chr_start
  | `Antisense -> r.chr_end
  in 
  r.chr, make pos pos
)

let of_file fn = 
  (File.lines_of fn |> Enum.skip 1)
  /@ Tsv.row_of_string
  /@ parse_row




















