open Batteries
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
   

let range_pos rx ry = min 0 (Range.gap rx ry)

let position map ((_,rx) as loc_x) = 
  let ((_,ry) as loc_y), _, _ = LMap.closest loc_x map in 
  range_pos rx ry

let histogram map xs = 
  Enum.fold
    (fun hist ((_,rx) as loc_x) -> 
      try 
        let pos = position map loc_x in 
        Histogram.increment hist pos 
      with Not_found -> hist)
    (Histogram.make compare histogram_bounds)
    xs

let ( !? ) x f = Option.bind x f

let histogram_count hist x = Histogram.(
  (find_bin_index hist x) !? (count hist)
)
               

let score chrom_size fa a fb b =
  let a = List.of_enum a
  and b = List.of_enum b in
  let amap = LMap.of_enum (List.enum a /@ (fun x -> fa x, x)) in 
  let hist_real = histogram amap (List.enum b /@ fb)
  and hist_rand = histogram amap (List.enum b /@ fb) in
  List.enum a
  /@ (fun (chra,ra) ->
    List.enum b
    /@ (fun (chrb,rb) -> 
      if chra <> chrb then Some 0.
      else (
        let pos = range_pos rx ry in
        match (histogram_count hist_real pos,
               histogram_count hist_rand pos) with
        | Some nreal, Some nrand -> 
            let r = (nreal -. nrand) /. nreal in 
            if classify_float r = FP_normal then Some r
            else None
        | _ -> None
      ))
  assert false




















