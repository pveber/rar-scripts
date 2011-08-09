open Batteries
open Oregon

let archive = Misc.unvac 
  ~host:"titan.ibisc.fr" 
  ~adr:"vac:bec7ce44a0702d261c56f88e31ec2dd6a1128f0d"
  (`directory ())

(* Fichier des pics sélectionnés *)
(* type col selected_peak = { *)
(*   chr : string ; *)
(*   st : int ; *)
(*   ed : int *)
(* } *)
(* let selected_peaks =  *)
(*   in_directory  *)
(*     archive ("clustering" et plus)  *)
(*     Tsv.({ has_header = assert false ; *)
(* 	   parse = Selected_peak.of_array }) *)

(*
let ( !! ) = Lazy.force

let nbclusters = lazy (
  shell <:sprint<ls $s:archive#path$/clustering/$s:assert false$/*.??? | wc -l>>
  |> int_of_string
)
*)

(* let clusters archive =  *)
(*   (in_directory archive ("clustering" et plus) ())#path *)
(*   |> File.lines_of  *)
(*   |> Enum.map int_of_string  *)
(*   |> Array.of_enum *)

(* let count_nb_clusters c =  *)
(*   Array.enum c *)
(*   |> PMap.of_enum *)
(*   |> PMap.cardinal *)

(* let cluster_homogeneity_pdf = F.make  *)
(*   (object *)
(*      method id = "Justin_chipseq_clustering" *)
(*      method deps = [] ++ archive *)
(*      method build path =  *)
(*        let rp = R.make ()  *)
(*        and clusters = clusters archive in *)
(*        let n = count_nb_clusters clusters in *)
(*        R.pdf rp path ; *)
(* 	 // préparer le plot en n fenêtres *)
(*        for i = 1 to n do *)
	 
(*        done *)
(*        R.devoff rp ; *)
(*        R.close rp *)
(*      method ty = `pdf *)
(*    end) *)

(* let filter_by_cluster a i =  *)
(*   let clusters = clusters () in *)
(*   Array.filteri (fun  *)
