  (* IDENTIFIER DES LOCUS SUPPORTANT LES PICS A TRAVERS DIFFERENTES CONDITIONS *)
open Batteries
open Printf
open Oregon
open Genome
open Oregon.Target
open Oregon.Target.Infix

module HRel = Accu.HRel

(*
let peaks_by_chromosome peaks = 
  HRel.of_enum
    (fun p -> p.PeakAnnotation.loc.Location.chr, p)
    (Array.enum peaks)

let summit p = PeakAnnotation.(p.loc.Location.st + p.summit)

let peak_islands peaks = PeakAnnotation.(
  let peaks = Array.of_list peaks in
  Array.sort (fun p q -> compare (summit p) (summit q)) peaks ;
  let current = ref [ peaks.(0) ] and r = ref [] in
  for i = 1 to Array.length peaks - 1 do 
    if summit peaks.(i) - summit (List.hd !current) > 100 then (
      r := !current :: !r ;
      current := [ peaks.(i) ]
    )
    else 
      current := peaks.(i) :: !current
  done ;
  r := !current :: !r ;
  List.rev_map List.rev !r
)

*)


let list_mean l = (List.fold_left ( + ) 0 l) / (List.length l)
let list_min l = List.fold_left min max_int l
let list_max l = List.fold_left max min_int l

(*
let peak_average l = 
  let summits = List.map summit l in
  let p = list_mean summits in
  Location.make 
    (List.hd l).PeakAnnotation.loc.Location.chr
    p p,
  list_max summits - list_min summits,
  List.length summits
*)    

let region chr island = 
  (chr,
   list_mean island,
   list_max island - list_min island,
   List.length island)

let islands chr summits = 
  let summits = Array.of_list summits in
  Array.sort compare summits ;
  let current = ref [ summits.(0) ] and r = ref [] in
  for i = 1 to Array.length summits - 1 do 
    if summits.(i) - (List.hd !current) > 100 then (
      r := !current :: !r ;
      current := [ summits.(i) ]
    )
    else 
      current := summits.(i) :: !current
  done ;
  r := !current :: !r ;
  List.rev_map (region chr) !r

let make ~summit ~location peaks = 
  let peaks_by_chr : (string, int) HRel.t = 
    HRel.of_enum 
      (fun p -> 
	 let l_p = location p in 
	 Location.(l_p.chr, l_p.st + summit p)) 
      (List.enum peaks) in
  Enum.fold 
    (fun accu (chr,summits) -> (islands chr summits) @ accu) 
    [] (HRel.enum peaks_by_chr)



let of_macs_targets ~radius ?(pvalue = 30.) macs_outputs = V.make 
  (object
     method id = sprintf "Binding_loci.of_macs_targets[r1,radius=%d,pvalue=%f]" radius pvalue
     method deps = [] ++* macs_outputs
     method build = 
       let all_peaks =
         ((List.enum macs_outputs)
             /@ Tsv.enum
             |> Enum.concat)
	 // Macs.(fun p -> p#pvalue > pvalue)
	 |> List.of_enum
       in
       make ~summit:(fun p -> p#summit) ~location:(fun l -> l#loc) all_peaks
       |> List.map (fun (chr,i,_,_) -> Location.make chr (i - radius) (i + radius))
   end)

  




  
