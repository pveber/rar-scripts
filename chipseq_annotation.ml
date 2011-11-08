open Batteries
open Oregon
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

type 'a annot = {
  coverage : ('a, int) HFun.t ;
  rpkm : ('a, float) HFun.t ;
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
    method id = "Chipseq_annotation[r1]" 
    method deps = [] ++ locations ++* (List.map chrmap all_samples) ++* (List.map nbmappings all_samples)
    method build = 
      let locations = locations#value in 
      let _ = Pmt.prime_cache() in
      let read_maps = HFun.cache (fun x -> (chrmap x)#value) in
      let library_size = HFun.cache (fun x -> (nbmappings x)#value) in
      let coverage x loc = Chr_map.nbregions_in (read_maps x) loc in
      let rpkm x loc = 
	float (coverage x loc) /. float (library_size x) *. 1e6 /. (loc |> Location.length |> float) *. 1000.
      and pvalue x y loc =
	pmt_test (coverage x loc) (coverage y loc) (library_size x) (library_size y) in
      let hfun f = 
	let init = List.enum design /@ fst 
	and f x = f x (List.assoc x design) in
	HFun.(make ~init f |> close) in
      Array.map 
	(fun loc -> { coverage =  hfun (fun x _ -> coverage x loc) ;
		      rpkm = hfun (fun x _ -> rpkm x loc) ;
		      pvalue = hfun (fun x y -> pvalue x y loc) })
	locations
  end)
