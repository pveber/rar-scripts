open Batteries
open Guizmin

type output

let closest_gene_for_each_bregion () = Resources.(
  (panRAR_regions ()
      |> Region_assoc.closest identity (Backup.get tss_map) )
  /@ (fun ((_,r_region) as region, ((_,r_tss),(_,_,dir)) as transcript) -> 
    let (_,r_region_center) = Location.center region in
    let pos = Region_assoc.stranded_range_pos ~from:(r_tss, dir) r_region_center in
    (region, transcript, pos))
)

let distance_breaks = [ -100000;-50000;-20000;-10000;-5000;-2000;-1000;0;1000;2000;5000;10000;20000;50000;100000]
let r_distance_breaks = 
  `l (R.floats (List.map Pervasives.float distance_breaks))

let rel_pos_bregions_fig_path = 
  "relative_position_of_bindings_wrt_closest_tss.pdf"

let rel_pos_bregions_fig_script = "\
rel_pos_bregions_fig <- function(path,
"

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

let app = Guizmin.d0 ("rar.app_region_assoc[r5]", []) (fun path ->
  ignore (Sys.command ("mkdir -p " ^ path)) ;
  (*  *)
  let closest_gene_for_each_bregion_list = List.of_enum (closest_gene_for_each_bregion ()) in
  let positions = List.map (Tuple3.third |- float_of_int) closest_gene_for_each_bregion_list in
  job path positions ;
  debug_bck path positions 

)

let rel_pos_bregions_fig app =
  select app rel_pos_bregions_fig_path




















