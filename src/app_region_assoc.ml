open Guizmin

type output
(*
let closest_gene_for_each_bregion () = Resource.(
  panRAR_regions ()
  |> Region_assoc.closest identity (Backup.get tss_map)
  /@ (fun ((_,r_region) as region, ((_,r_tss),_) as transcript)  -> 
    let pos = Region_assoc.range_pos r_region r_tss)

)

let make_
*)

let app = Guizmin.d0 ("rar.app_region_assoc[r1]", []) (fun path ->
  assert false
)

let relative_position_of_bindings_wrt_gene_tss_fig app =
  select app "relative_position_of_bindings_wrt_gene_tss.pdf"




















