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

let app = Guizmin.d0 ("rar.app_region_assoc[r1]", []) (fun path ->
  ignore (Sys.command ("mkdir -p " ^ path)) ;
  (*  *)
  let closest_gene_for_each_bregion_list = List.of_enum (closest_gene_for_each_bregion ()) in
  RgrDevices.pdf (Filename.concat path rel_pos_bregions_fig_path) ;
  Rgraphics.hist 
    ~breaks:r_distance_breaks
    ~main:"Relative position of binding regions wrt to closest TSS"
    ~xlab:"Relative position (in bp)"
    ~ylab:"Frequence"
    (closest_gene_for_each_bregion_list
       |> List.map (Tuple3.third |- float_of_int)
       |> List.filter (fun x -> abs_float x < 100000.)
       |> R.floats) |> ignore ;
  RgrDevices.dev_off () ;

)

let rel_pos_bregions_fig app =
  select app rel_pos_bregions_fig_path




















