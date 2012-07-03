(* let () = *)
(*   let gene_and_neighbours = List.of_enum (gene_neighbouring_regions ()) in   *)
(*   let open Rnaseq_table in *)
(*       List.iter (fun (gene, regions) -> *)
(*         printf "%s %s %d %d %s\n%!" gene.gene_symbol gene.chr gene.chr_start gene.chr_end (Location.to_string (tss_loc_of_row gene)) ; *)
(*         Array.iter *)
(*           (fun (loc,_) -> print_endline (Location.to_string loc)) *)
(*           regions *)
(*       ) *)
(*         gene_and_neighbours *)

(*
let () = 
  gene_test_modulated_vs_region_in_promoter ()
*)


(* let () =  *)
(*   ignore (Backup.get closest_gene_from_region) *)

(* let () = *)
(*   gene_test_modulated_vs_region_in_promoter () *)

(* let () = *)
(*   gene_test_modulated_vs_motif (`direct, 0) *)



(* let () = *)
(*   ignore (Backup.get closest_gene_from_region) *)

let _ = Guizmin.eval App_region_assoc.(relative_position_of_bindings_wrt_gene_tss_fig app)










