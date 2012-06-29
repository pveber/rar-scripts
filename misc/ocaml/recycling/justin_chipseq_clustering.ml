open Batteries
open Printf
open Oregon
open Target.Infix

(* sortie par justin *)
type col selected_peak = {
  chr : string ;
  st : int ;
  ed : int ;
  t0 : float ;
  t2 : float ;
  t24 : float ;
  t48 : float ;
  cl2 : int ;
  cl3 : int ;
  cl4 : int ;
  cl5 : int ;
  cl6 : int ;
  cl7 : int ;
  cl8 : int ;
  cl9 : int ;
  cl10 : int 
}

let s_of_tf = function `rar -> "rar" | `rxr -> "rxr"
let s_of_cm = function `spectral -> "spectral" | `hsic -> "hsic"

let get_cluster n p = 
  try
  (Selected_peak.to_array p).(Selected_peak.Index.cl2 + n - 2)
  |> int_of_string
  with _ -> failwith "smartass"

let clustering_path tf cm = 
  sprintf "extdata/justin/clustering/chipseq/%s-%s.txt" (s_of_tf tf) (s_of_cm cm)

let load_peaks tf cm = 
  File.lines_of (clustering_path tf cm)
  |> Enum.map (Sle.SleString.split '\t')
  |> Enum.map Selected_peak.of_array
  |> Array.of_enum

let count_nb_clusters c =
  Array.enum c
  |> PMap.of_enum
  |> PMap.enum 
  |> Enum.hard_count

let cluster_homogeneity_fig tf cm n path = 
  let rp = R.make () 
  and peaks = load_peaks tf cm in
  let cluster_assignment = Array.map (get_cluster n) peaks in
  R.pdf rp path ;
  for j = 0 to n - 1 do
    R.c rp "plot(c(0,48),c(-5,5),type='n')" ;
    Sle.SleArray.filteri (fun i _ -> cluster_assignment.(i) = j) peaks 
    |> Array.iter (fun p -> Selected_peak.(R.c rp "lines(c(0,2,24,48),c(%f,%f,%f,%f),lwd=0.3)" p.t0 p.t2 p.t24 p.t48)) ;
  done ;
  R.devoff rp ;
  R.close rp

(* let cluster_homogeneity_fig_target tf cm n = Target.F.make *)
(*   (object *)
(*      method id = "Justin_chipseq_clustering.cluster_homogeneity_fig[r2]" *)
(*      method deps = [] ++ archive *)
(*      method build path = cluster_homogeneity_fig tf cm n path *)
(*      method ty = `pdf *)
(*    end) *)

(* let filter_by_cluster a i =  *)
(*   let clusters = clusters () in *)
(*   Array.filteri (fun  *)
