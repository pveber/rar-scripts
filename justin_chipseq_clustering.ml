open Batteries
open Printf
open Oregon
open Target.Infix

let archive = Misc.unvac 
  ~host:"titan.ibisc.fr" 
  ~adr:"vac:bec7ce44a0702d261c56f88e31ec2dd6a1128f0d"
  (`directory ())

let prepare_archive () = 
  Sle.sh "cd %s && find . -name '*.gz' -exec gunzip {} \\;" archive#path

(* Fichier des pics sélectionnés *)
type col selected_peak = {
  chr : string ;
  st : int ;
  ed : int ;
  t0 : float ;
  t2 : float ;
  t24 : float ;
  t48 : float
}

let selected_peaks =
  Target.in_directory
    archive "clustering/res/chipseq/data.txt"
    Tsv.({ has_header = true ;
	   parse = Selected_peak.of_array })



(*
let ( !! ) = Lazy.force

let nbclusters = lazy (
  shell <:sprint<ls $s:archive#path$/clustering/$s:assert false$/*.??? | wc -l>>
  |> int_of_string
)
*)
    
let cluster_assignment tf n =
  let tf = match tf with `rar -> "rar" | `rxr -> "rxr" in
  (Target.in_directory archive (sprintf "clustering/res/chipseq/hsic/%s/%s-%d" tf tf n) ())#path
  |> File.lines_of
  |> Enum.skip 1
  |> Enum.map int_of_string
  |> Array.of_enum

let count_nb_clusters c =
  Array.enum c
  |> PMap.of_enum
  |> PMap.enum 
  |> Enum.hard_count

let cluster_homogeneity_fig tf n = Target.F.make
  (object
     method id = "Justin_chipseq_clustering.cluster_homogeneity_fig[r1"
     method deps = [] ++ archive
     method build path =
       let rp = R.make ()
       and cluster_assignment = cluster_assignment tf n 
       and peaks = Tsv.enum selected_peaks |> Array.of_enum in
       let sn = int_of_float (ceil (sqrt (float n))) in
       R.pdf rp path ;
       R.c rp "par(mfrow=c(%d,%d))" sn sn ;
       for j = 1 to n do
	 R.c rp "plot(c(1,3),c(-5,5),type='n')" ;
	 Sle.SleArray.filteri (fun i _ -> cluster_assignment.(i) = j) peaks 
	 |> Array.iter (fun p -> Selected_peak.(R.c rp "lines(c(1,2,3),c(%f,%f,%f),lwd=0.3)" p.t2 p.t24 p.t48)) ;
       done ;
       R.devoff rp ;
       R.close rp
     method ty = `pdf
   end)

(* let filter_by_cluster a i =  *)
(*   let clusters = clusters () in *)
(*   Array.filteri (fun  *)
