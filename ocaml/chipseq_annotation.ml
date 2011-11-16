open Batteries
open Oregon
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

type 'a annot = {
  loc : Location.t ;
  counts : ('a, int) HFun.t ;
  coverage : ('a, float) HFun.t ;
  pvalue : ('a, float) HFun.t
}

let pmt_test x1 x2 r1 r2 = 
  let l1 = float x1 /. float r1 
  and l2 = float x2 /. float r2 in 
  if l1 < l2 then Pmt.log_poisson_margin x1 x2 r1 r2
  else Pmt.log_poisson_margin x2 x1 r2 r1

let make design locations = 
  let all_samples = 
    List.fold_left (fun accu (x,y) -> x :: y :: accu) [] design
    |> List.sort_unique compare 
  and bowtie = HFun.cache B.Chipseq.bowtie_wodup in
  let chrmap = HFun.cache (bowtie |- Sam.to_chr_map)
  and nbmappings = HFun.cache (bowtie |- Sam.nbmappings) in
  V.make (object
    method id = "Chipseq_annotation[r2]" 
    method deps = [] ++ locations ++* (List.map chrmap all_samples) ++* (List.map nbmappings all_samples)
    method build = 
      let locations = locations#value in 
      let _ = Pmt.prime_cache() in
      let read_maps = HFun.cache (fun x -> (chrmap x)#value) in
      let library_size = HFun.cache (fun x -> (nbmappings x)#value) in
      let counts x loc = Chr_map.nbregions_in (read_maps x) loc in
      let coverage x loc = 
	float (counts x loc) /. float (library_size x) *. 1e6 /. (loc |> Location.length |> float) *. 1000.
      and pvalue x y loc =
	-. pmt_test (counts x loc) (counts y loc) (library_size x) (library_size y) in
      let hfun f = 
	let init = List.enum design /@ fst 
	and f x = f x (List.assoc x design) in
	HFun.(make ~init f |> close) in
      Array.map 
	(fun loc -> { loc      = loc ;
		      counts   =  hfun (fun x _ -> counts x loc) ;
		      coverage = hfun (fun x _ -> coverage x loc) ;
		      pvalue   = hfun (fun x y -> pvalue x y loc) })
	locations
  end)

let domain annot = match annot with 
    [| |] -> []
  | t -> HFun.domain t.(0).counts |> List.of_enum
  
let tsv output_path annot = 
  let domain = domain annot in
  let fields f hfun = 
    (List.enum domain) /@ (( $ ) hfun |- f) |> Array.of_enum
  and labels l = 
    (List.enum domain) /@ (fun x -> l ^ " " ^ (B.Chipseq.string_of_sample x)) |> Array.of_enum
  in
  Array.enum annot
  /@ Location.(fun a -> Array.concat [
    [| a.loc.chr ; string_of_int a.loc.st ; string_of_int a.loc.ed |] ;
    fields string_of_int a.counts ;
    fields string_of_float a.coverage ;
    fields string_of_float a.pvalue 
  ])
  |> Enum.(append (singleton (
       Array.concat [
	 [| "chrom" ; "start" ; "end" |] ;
	 labels "counts" ;
	 labels "coverage" ;
	 labels "pvalue" 
       ])))
  |> Enum.map Tsv.string_of_line
  |> File.write_lines output_path

