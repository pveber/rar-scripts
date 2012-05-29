open Batteries
open Printf
open Biocaml
open GenomeMap

let closest f map regions = 
  let aux x = 
    let loc_y, y, _ = LMap.closest (f x) map in 
    x, (loc_y,y)
  in 
  Enum.map aux regions

let histogram_bounds = [ 
  -100000; -50000; -20000; -10000 ; -5000 ;
  -2000 ; -1000 ; 0 ; 1000; 2000; 5000;
  10000; 20000; 50000; 100000 
]
   
let range_overlap s s' = Range.(
  let p s s' = (s.hi >= s'.lo) && (s.hi <= s'.hi) in 
  p s s' || p s' s
)

let range_pos ~from r = Range.(
  if range_overlap from r then 0
  else
    let a, b = r.hi - from.lo, r.lo - from.hi in
    if abs a < abs b then a else b
)

let position map ((_,rx) as loc_x) = 
  let (_,ry), _, _ = LMap.closest loc_x map in 
  range_pos rx ry

let histogram map xs = 
  Enum.fold
    (fun hist ((chr,rx) as loc_x) -> 
      try 
        let pos = position map loc_x in 
        Histogram.increment hist pos
      with Not_found -> hist)
    (Histogram.make compare histogram_bounds)
    xs

let ( |?> ) x f = Option.bind f x

let histogram_count hist x = Histogram.(
  Option.map
    (count hist)
    (find_bin_index hist x)
)
               
let print_histogram hist = 
  List.iter
    (fun ((a,b),c) -> printf "%d %d %f\n%!" a b c)
    (Histogram.to_list hist)

let score chrom_size fa a fb b =
  let a = List.of_enum a
  and b = List.of_enum b in
  let amap = LMap.of_enum (List.enum a /@ (fun x -> fa x, x)) in 
  let hist_real = histogram amap (List.enum b /@ fb)
  and hist_rand = histogram amap (List.enum b /@ fb) in
  print_histogram hist_real ;
  print_histogram hist_rand ;
  List.enum a
  /@ (fun a ->
    let (chra,ra) = fa a in
    List.enum b
    //@ (fun b -> 
      let (chrb,rb) = fb b in
      if chra <> chrb then None
      else (
        let pos = range_pos ra rb in
        match (histogram_count hist_real pos,
               histogram_count hist_rand pos) with
        | Some nreal, Some nrand -> 
            let r = (nreal -. nrand) /. nreal in 
            (* printf "%f\n%!" r ; *)
            if classify_float r = FP_normal && r > 0. 
            then Some (a,b,r)
            else None
        | _ -> None
      )))
  |> Enum.concat




















