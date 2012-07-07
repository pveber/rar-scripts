open Batteries
open Guizmin
open Guizmin_bioinfo.MBSchema

type output

let closest_gene_for_each_bregion () = Resources.(
  (panRAR_regions ()
      |> Region_assoc.closest identity (Backup.get tss_map) )
  /@ (fun (region, (_,transcript)) -> 
    let region_center = Location.center region in
    let pos = Transcript.position2tss transcript region_center in
    (region, transcript, pos))
)


(** for each gene, find a binding region whose to distance to one of
    the gene tss is minimal *)
let closest_bregion_for_each_gene () =

  (* an accumulator that keeps the closest transcript seen so far for each gene *)
  let accu = Biocaml_accu.create 
    None identity
    (fun (transcript, (bregion,bregion_id)) current -> 
      let current_is_closer = match current with
        | None -> false 
        | Some (transcript',_,_) -> 
            let d  = Location.dist (Transcript.tss transcript)  bregion
            and d' = Location.dist (Transcript.tss transcript') bregion in
            d' < d
      in
      if current_is_closer then current
      else Some (transcript,bregion,bregion_id))
  in 
  let bregions_map = Resources.panRAR_region_summits_map () in 
  let transcript2closest_bregion = 
    Region_assoc.closest Transcript.tss bregions_map (Ensembl.transcripts_enum (Resources.gtf ())) 
  in
  Enum.iter 
    (fun ((transcript, _) as x) -> Biocaml_accu.add accu transcript.Transcript.gene_id x)
    transcript2closest_bregion ;
  let hastbl = Hashtbl.of_enum (Biocaml_accu.enum accu) in
  
  assert false


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
  select app "relative_position_of_bindings_wrt_closest_tss.pdf"




















