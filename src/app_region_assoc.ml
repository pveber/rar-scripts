open Batteries
open Biocaml.GenomeMap
open Guizmin
open Guizmin_bioinfo.MBSchema

(* hack to load libraries in the R interpreter *)
let _ = Rgraphics.hist, RgrDevices.dev_off

type output

let fun_of_enum f g e = 
  let tbl = Hashtbl.of_enum (e /@ (fun x -> f x, g x)) in 
  fun x -> try Some (Hashtbl.find tbl x) with Not_found -> None

let fun_of_enum_exn f g e = 
  let tbl = Hashtbl.of_enum (e /@ (fun x -> f x, g x)) in 
  fun x -> Hashtbl.find tbl x

let closest_gene_for_each_bregion () = Resources.(
  (panRAR_regions ()
      |> Region_assoc.closest identity (Backup.get tss_map) )
  /@ (fun (region, (_,transcript)) -> 
    let region_center = Location.center region in
    let pos = Transcript.position2tss transcript region_center in
    (object
      method region = region
      method transcript = transcript
      method relpos2transcript = pos
     end))
)


(** for each gene, find a binding region whose to distance to one of
    the gene tss is minimal *)
let closest_bregion_for_each_gene () =
  let bregions_map = Resources.panRAR_region_summits_map () in 
  let foreach_transcript t =
    let loc, idx, d = LMap.closest (Transcript.tss t) bregions_map in 
    d, loc, idx, t in
  let foreach g =
    try 
      let _, loc, idx, t = 
        List.map foreach_transcript g.Gene.transcripts
        |> List.reduce min
      in
      Some (object
        method gene = g 
        method transcript = t 
        method region_loc = loc
        method region_idx = idx
        method position2tss = Transcript.position2tss t loc
       end)
    with Not_found -> None
  in 
  Ensembl.genes_enum (Resources.gtf ())
  //@ foreach

let job_script = "\
job <- function(path, closest_tss_rel_pos) {
  breaks <- c(-1000000,-100000,-50000,-20000,-10000,-5000,-2000,-1000,0,1000,2000,5000,10000,20000,50000,100000,1000000)
  breaks_ranges <- abs(breaks[1:length(breaks) - 1] - breaks[2:length(breaks)])
  breaks_legend <- c('[-1000000;-100000]','[-100000;-50000]','[-50000;-20000]','[-20000;-10000]','[-10000;-5000]','[-5000;-2000]','[-2000;-1000]','[-1000;0]','[0;1000]','[1000;2000]','[2000;5000]','[5000;10000]','[10000;20000]','[20000;50000]','[50000;100000]','[100000;1000000]')
  rel_pos_hist <- hist(closest_tss_rel_pos[abs(closest_tss_rel_pos) < 1000000], breaks,plot=F)
  pdf(paste(path,'relative_position_of_bindings_wrt_closest_tss.pdf',sep='/'))
  barplot(rel_pos_hist$counts / breaks_ranges, names.arg=breaks_legend,  las=2, cex.names=0.5)
  dev.off()

}
"
let job = 
  let _ = R.eval_string job_script in 
  let symbol = R.symbol "job" in
  fun path closest_tss_rel_pos ->
    R.eval symbol [
      R.arg R.string path ;
      R.arg R.floats closest_tss_rel_pos ;
     ]
  |> ignore

let debug_bck_script = "\
debug_bck <- function(path, closest_tss_rel_pos) {
  backup <- list(closest_tss_rel_pos = closest_tss_rel_pos)
  save(backup,file=paste(path,'debug.Rd',sep='/'))
}
"
let debug_bck = 
  let _ = R.eval_string debug_bck_script in 
  let symbol = R.symbol "debug_bck" in
  fun path closest_tss_rel_pos ->
    R.eval symbol [
      R.arg R.string path ;
      R.arg R.floats closest_tss_rel_pos ;
    ]
  |> ignore


(** TSV associant pour chaque région le gène le plus proche *)
let closest_gene_from_region_tsv path closest_gene_from_region_list =
  List.enum closest_gene_from_region_list
  /@ (fun r -> [| 
    Location.chr r#region ; 
    string_of_int (Location.st r#region) ;
    string_of_int (Location.ed r#region) ;
    r#transcript.Transcript.gene_id ;
    r#transcript.Transcript.id ;
    string_of_int r#relpos2transcript ;
  |])
  /@ Tsv.string_of_row
  |> File.write_lines (Filename.concat path "closest_gene_from_region.tsv")

(** Figure montrant la proportion de chaque classe d'expression en
    fonction de la distance à la plus proche région *)
let fig_expr_class_dist_wrt_closest_bregion_pos 
    path genes closest_bregion_of_gene_id  expression_of_gene_id
    upregulated downregulated expressed time  = 
  let expression_class r =
    if upregulated r then 3
    else if downregulated r then 2
    else if expressed r then 1
    else 0
  in 
  let positions, expr_class = 
    (List.enum genes)
    //@ (fun g -> 
      try Some (closest_bregion_of_gene_id g.Gene.id,
                expression_class (expression_of_gene_id g.Gene.id))
      with Not_found -> None)
    |> List.of_enum
    |> List.split
  in 
  let script = Printf.sprintf "\
f <- function(path, position, expr_class) {
  backup <- list(position = position, expr_class = expr_class)
  save(backup,file=paste(path,'expr_class_dist_wrt_closest_bregion_pos.Rd',sep='/'))
  breaks <- c(-1000000,-100000,-50000,-20000,-10000,-5000,-2000,-1000,0,1000,2000,5000,10000,20000,50000,100000,1000000)
  breaks_ranges <- abs(breaks[1:length(breaks) - 1] - breaks[2:length(breaks)])
  x <- table(expr_class,cut(position,breaks))
  pdf(paste(path,'relative_position_closest_bregion.pdf',sep='/'))
  barplot(table(cut(position,breaks)) / breaks_ranges, las=2, cex.names=0.5)
  dev.off()
  pdf(paste(path,'expr_class_dist_wrt_closest_bregion_pos%s.pdf',sep='/'))
  barplot(t(t(x) / colSums(x)), las=2, cex.names=0.5)
  dev.off()
}
" time 
  in
  let _ = R.eval_string script in 
  R.eval (R.symbol "f") [
    R.arg R.string path ;
    R.arg R.ints positions ;
    R.arg R.ints expr_class
  ]
  |> ignore
  

let fig_expr_class_prop_wrt_closest_bregion_pos path genes closest_bregion_of_gene_id expression_of_gene_id = 
  let open Rnaseq_table in
  let p x = Filename.concat path (x ^ "_prop_wrt_closest_bregion_pos.pdf") in
  let available_ids, positions = 
    (List.enum genes)
    //@ (fun g -> 
      try Some (g.Gene.id, closest_bregion_of_gene_id g.Gene.id)
      with Not_found -> None)
    |> List.of_enum
    |> List.split
  in
  let vec p = List.map (expression_of_gene_id |- p)  available_ids in
  let script = Printf.sprintf "\
f <- function(path, position, freq_6H,freq_12H,freq_24H,freq_36H,freq_48H) {
  backup <- list(position = position, freq_6H = freq_6H, freq_12H = freq_12H, freq_24H = freq_24H, freq_36H = freq_36H, freq_48H = freq_48H)
  save(backup,file=paste(path,'.Rd',sep=''))
  breaks <- c(-1000000,-100000,-50000,-20000,-10000,-5000,-2000,-1000,0,1000,2000,5000,10000,20000,50000,100000,1000000)
  breaks_ranges <- abs(breaks[1:length(breaks) - 1] - breaks[2:length(breaks)])
  breaks_legend <- c('[-1000000;-100000]','[-100000;-50000]','[-50000;-20000]','[-20000;-10000]','[-10000;-5000]','[-5000;-2000]','[-2000;-1000]','[-1000;0]','[0;1000]','[1000;2000]','[2000;5000]','[5000;10000]','[10000;20000]','[20000;50000]','[50000;100000]','[100000;1000000]')
  tt_6H <- table(cut(position,breaks),freq_6H)
  xx_6H <- tt_6H / rowSums(tt_6H)
  tt_12H <- table(cut(position,breaks),freq_12H)
  xx_12H <- tt_12H / rowSums(tt_12H)
  tt_24H <- table(cut(position,breaks),freq_24H)
  xx_24H <- tt_24H / rowSums(tt_24H)
  tt_36H <- table(cut(position,breaks),freq_36H)
  xx_36H <- tt_36H / rowSums(tt_36H)
  tt_48H <- table(cut(position,breaks),freq_48H)
  xx_48H <- tt_48H / rowSums(tt_48H)
  ylim <- c(0,1.1 * max(rbind(xx_6H[,2],xx_12H[,2],xx_24H[,2],xx_36H[,2],xx_48H[,2])))
  x <- 1:length(breaks_legend)
  pdf(path)
  plot(x,xx_6H[,2],type='l',ylim=ylim,xaxt='n',col=rainbow(5)[1],xlab='',ylab='')
  axis(1,at=1:length(breaks_legend),labels=breaks_legend,las=2,cex.axis=0.5)
  lines(x,xx_12H[,2],type='l',col=rainbow(5)[2])
  lines(x,xx_24H[,2],type='l',col=rainbow(5)[3])
  lines(x,xx_36H[,2],type='l',col=rainbow(5)[4])
  lines(x,xx_48H[,2],type='l',col=rainbow(5)[5])
  legend('topright',legend=c('6H','12H','24H','36H','48H'),fill=rainbow(5))
  dev.off()
}
" 
  in
  let _ = R.eval_string script in 
  let f_stub = (R.symbol "f") in
  let script path freq_6H freq_12H freq_24H freq_36H freq_48H = 
    R.eval f_stub [
      R.arg R.string path ;
      R.arg R.ints positions ;
      R.arg R.bools freq_6H ;
      R.arg R.bools freq_12H ;
      R.arg R.bools freq_24H ;
      R.arg R.bools freq_36H ;
      R.arg R.bools freq_48H ;
    ]
  in 
  ignore (script (p "upregulated") (vec upregulated_6H) (vec upregulated_12H) (vec upregulated_24H) (vec upregulated_36H) (vec upregulated_48H)) ;
  ignore (script (p "downregulated") (vec downregulated_6H) (vec downregulated_12H) (vec downregulated_24H) (vec downregulated_36H) (vec downregulated_48H))


let comodulation_under_common_bregion expr_pred1 expr_pred2 genes expression_of_gene_id =
  let radius = 20000 in
  let region_map = Resources.panRAR_regions_map () in
  let tss_of_gene g = List.map Transcript.tss g.Gene.transcripts in
  let compare_genes g1 g2 = compare g1.Gene.id g2.Gene.id in
  let gene_promoter g = 
    List.enum (tss_of_gene g) 
    /@ (Location.relmove (-radius) radius)
    |> Selection.of_enum in
  let genes_common_promoter g g' =
    Selection.inter (gene_promoter g) (gene_promoter g')
  in
  let common_region g g' = 
    genes_common_promoter g g'
    |> Selection.enum 
    |> Enum.exists
        (fun loc -> LMap.intersects loc region_map)
  in
  let classify (g1,g2) =
    let pred1 g = expr_pred1 (expression_of_gene_id g.Gene.id)
    and pred2 g = expr_pred2 (expression_of_gene_id g.Gene.id) in
    ((pred1 g1 && pred2 g2) || (pred1 g2 && pred2 g1)),
    common_region g1 g2
  in
  Region_assoc.gene_tss_proximity_graph compare_genes tss_of_gene ~radius (List.enum genes)
  /@ (fun (g, neighbours) ->
    (List.enum neighbours)
    //@ (fun g' -> if compare_genes g g' < 0 then Some (g,g') else None))
  |> Enum.concat
  |> Biocaml_accu.counts classify 
  |> Enum.iter (fun ((b1, b2), c) -> Printf.printf "%b\t%b\t%d\n" b1 b2 c)


let app = Guizmin.d0 ("rar.app_region_assoc[r19]", []) (fun path ->
  ignore (Sys.command ("mkdir -p " ^ path)) ;
  (*  *)
  let closest_gene_for_each_bregion_list = List.of_enum (closest_gene_for_each_bregion ()) in
  let closest_gene_positions = 
    List.map 
      (fun r -> float_of_int r#relpos2transcript) 
      closest_gene_for_each_bregion_list in
  let genes = List.of_enum (Ensembl.genes_enum (Resources.gtf ())) in
  let expression_of_gene_id = 
    fun_of_enum_exn
      (fun r -> r.Rnaseq_table.ensembl_gene_id)
      identity
      (Resources.rnaseq_table ())
  in 
  Rnaseq_table.(
    print_endline "Co-upregulation 6H" ;
    comodulation_under_common_bregion upregulated_6H upregulated_6H genes expression_of_gene_id ;
    print_endline "Co-upregulation 12H" ;
    comodulation_under_common_bregion upregulated_12H upregulated_12H genes expression_of_gene_id ;
    print_endline "Co-upregulation 24H" ;
    comodulation_under_common_bregion upregulated_24H upregulated_24H genes expression_of_gene_id ;
    print_endline "Co-upregulation 48H" ;
    comodulation_under_common_bregion upregulated_48H upregulated_48H genes expression_of_gene_id ;
    print_endline "Co-downregulation 6H" ;
    comodulation_under_common_bregion downregulated_6H downregulated_6H genes expression_of_gene_id ;
    print_endline "Co-downregulation 12H" ;
    comodulation_under_common_bregion downregulated_12H downregulated_12H genes expression_of_gene_id ;
    print_endline "Co-downregulation 24H" ;
    comodulation_under_common_bregion downregulated_24H downregulated_24H genes expression_of_gene_id ;
    print_endline "Co-downregulation 48H" ;
    comodulation_under_common_bregion downregulated_48H downregulated_48H genes expression_of_gene_id ;
    print_endline "Co-updownregulation 6H" ;
    comodulation_under_common_bregion upregulated_6H downregulated_6H genes expression_of_gene_id ;
    print_endline "Co-updownregulation 12H" ;
    comodulation_under_common_bregion upregulated_12H downregulated_12H genes expression_of_gene_id ;
    print_endline "Co-updownregulation 24H" ;
    comodulation_under_common_bregion upregulated_24H downregulated_24H genes expression_of_gene_id ;
    print_endline "Co-updownregulation 48H" ;
    comodulation_under_common_bregion upregulated_48H downregulated_48H genes expression_of_gene_id ;
    exit 42 ;
    let closest_bregion_of_gene_id = 
      fun_of_enum_exn
        (fun r -> r#gene.Gene.id) 
        (fun r -> r#position2tss)
        (closest_bregion_for_each_gene ())
    in
    let f = fig_expr_class_dist_wrt_closest_bregion_pos path genes closest_bregion_of_gene_id expression_of_gene_id in 
    f upregulated downregulated expressed "" ;
    f upregulated_6H downregulated_6H expressed_6H "_6H" ;
    f upregulated_12H downregulated_12H expressed_12H "_12H" ;
    f upregulated_24H downregulated_24H expressed_24H "_24H" ;
    f upregulated_36H downregulated_36H expressed_36H "_36H" ;
    fig_expr_class_prop_wrt_closest_bregion_pos path genes closest_bregion_of_gene_id expression_of_gene_id 
  ) ;
  job path closest_gene_positions ;
  debug_bck path closest_gene_positions ;
  closest_gene_from_region_tsv path closest_gene_for_each_bregion_list ;
)

let fig_rel_pos_bregions app =
  select app "relative_position_of_bindings_wrt_closest_tss.pdf"

let fig_rel_pos_closest_bregion app =
  select app "relative_position_closest_bregion.pdf"

let fig_expr_class_dist_wrt_closest_bregion_pos app =
  select app "expr_class_dist_wrt_closest_bregion_pos.pdf"

let fig_expr_class_dist_wrt_closest_bregion_pos_6H app =
  select app "expr_class_dist_wrt_closest_bregion_pos_6H.pdf"

let fig_expr_class_dist_wrt_closest_bregion_pos_12H app =
  select app "expr_class_dist_wrt_closest_bregion_pos_12H.pdf"

let fig_expr_class_dist_wrt_closest_bregion_pos_24H app =
  select app "expr_class_dist_wrt_closest_bregion_pos_24H.pdf"

let fig_expr_class_dist_wrt_closest_bregion_pos_36H app =
  select app "expr_class_dist_wrt_closest_bregion_pos_36H.pdf"

let fig_upregulated_prop_wrt_closest_bregion_pos app = 
  select app "upregulated_prop_wrt_closest_bregion_pos.pdf"

let fig_downregulated_prop_wrt_closest_bregion_pos app = 
  select app "downregulated_prop_wrt_closest_bregion_pos.pdf"
  




















