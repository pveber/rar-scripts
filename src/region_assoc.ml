open Batteries
open Printf
open Biocaml
open GenomeMap
open Guizmin_bioinfo.MBSchema

let closest f map regions = 
  let aux x = 
    let loc_y, y, _ = LMap.closest (f x) map in 
    x, (loc_y,y)
  in 
  Enum.map aux regions

let neighbours f dmax map x =
  let loc_x = Location.relmove (-dmax) dmax (f x) in
  LMap.intersecting_elems loc_x map |> Array.of_enum 


let histogram_bounds = [ 
  -100000; -50000; -20000; -10000 ; -5000 ;
  -2000 ; -1000 ; 0 ; 1000; 2000; 5000;
  10000; 20000; 50000; 100000 
]
   

let position map loc_x = 
  let loc_y, _, _ = LMap.closest loc_x map in 
  Location.position ~from:loc_x loc_y

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
  let control = 
    let b = Array.of_list b in
    Random_regions.uniform ~chrom_size fb b (Array.length b) in 
  let hist_real = histogram amap (List.enum b /@ fb)
  and hist_rand = histogram amap control in
  print_histogram hist_real ;
  print_histogram hist_rand ;
  List.enum a
  /@ (fun a ->
    let loc_a = fa a in
    List.enum b
    //@ Location.(fun b -> 
      let loc_b = fb b in
      if chr loc_a  <> chr loc_b then None
      else (
        let pos = position loc_a loc_b in
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

let gene_re_graph tss_of_gene loc_of_re ~radius genes relts =
  let re_map = LMap.of_enum (relts /@ (fun x -> loc_of_re x, x)) in
  genes /@ (fun g ->
    let selected_relts = 
      List.map
        (fun tss_loc -> 
          let loc = Location.relmove (-radius) radius tss_loc in
          LMap.intersecting_elems loc re_map /@ snd |> Set.of_enum)
        (tss_of_gene g)
      |> List.fold_left Set.union Set.empty
    in 
    (g, selected_relts))


let gene_tss_proximity_graph (type gene) compare_genes tss_of_gene ~radius genes =
  let module G = struct
    type t = gene
    let compare = compare_genes
  end in
  let module GSet = Set.Make(G) in
  let genes = List.of_enum genes in
  let tss_map = 
    List.enum genes 
    /@ (fun g -> 
      let tss_g = tss_of_gene g in
      List.enum tss_g /@ (fun tss -> tss, g))
    |> Enum.concat
    |> LMap.of_enum
  in
  List.enum genes 
  /@ (fun g ->
    let neighbours = 
      (List.enum (tss_of_gene g))
      /@ (Location.relmove (-radius) radius)
      /@ (fun zone ->
        LMap.intersecting_elems zone tss_map /@ snd)
      |> Enum.concat
      |> GSet.of_enum
      |> GSet.elements
    in g, neighbours)

