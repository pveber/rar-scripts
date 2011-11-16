open Batteries
open Printf
open Sle.Infix
open Genome
open Oregon
open Target.Infix

open B.B

type t = {
  loc : Location.t ;
  coverage : (Chipseq.sample, int) Sle.hfun ;
  rpkm : (Chipseq.sample, float) Sle.hfun ;
  pvalue : (Chipseq.sample * Chipseq.sample, float) Sle.hfun ;
}

let samples = [ 
  `F9_WT_panRAR_1 ; 
  `F9_ATRA_panRAR_1 ;
  `F9_ATRA24_panRAR_1 ;
  `F9_ATRA48_panRAR_1 ;
  `F9_WT_panRXR_1 ; 
  `F9_ATRA2_panRXR_1 ;
  `F9_ATRA24_panRXR_1 ; 
  `F9_ATRA48_panRXR_1 ;
  `Input_1
]


let paired_samples = [
  `F9_WT_panRAR_1, `Input_1 ; 
  `F9_WT_panRXR_1, `Input_1 ; 
  `F9_ATRA_panRAR_1, `Input_1 ; 
  `F9_ATRA2_panRXR_1, `Input_1 ;
  `F9_ATRA24_panRAR_1, `Input_1 ; 
  `F9_ATRA24_panRXR_1, `Input_1 ;
  `F9_ATRA48_panRAR_1, `Input_1 ; 
  `F9_ATRA48_panRXR_1, `Input_1 ;
  `F9_WT_panRAR_1, `F9_ATRA_panRAR_1 ;
  `F9_WT_panRXR_1, `F9_ATRA2_panRXR_1 ;
  `F9_ATRA_panRAR_1, `F9_ATRA24_panRAR_1 ;
  `F9_ATRA2_panRXR_1, `F9_ATRA24_panRXR_1 ;
  `F9_ATRA24_panRAR_1, `F9_ATRA48_panRAR_1 ;
  `F9_ATRA24_panRXR_1, `F9_ATRA48_panRXR_1 ;
  `F9_WT_panRAR_1, `F9_ATRA24_panRAR_1 ;
  `F9_WT_panRXR_1, `F9_ATRA24_panRXR_1 ;
  `F9_WT_panRAR_1, `F9_ATRA48_panRAR_1 ;
  `F9_WT_panRXR_1, `F9_ATRA48_panRXR_1 ;
  `F9_ATRA_panRAR_1, `F9_ATRA48_panRAR_1 ;
  `F9_ATRA2_panRXR_1, `F9_ATRA48_panRXR_1 ;
]

module L = Labs.Make(B.B)
open B.B


let islands = L.binding_loci (List.filter ( (<>) `Input_1 ) samples)

let radius = 250

let location_of_island n i =
  let x,_,_ = L.peak_average i in 
  Location.(make x.chr (x.st - n) (x.ed + n))

module Annotation(X : sig end) = struct
  let _ = Pmt.prime_cache()
    
  let ( ! ) = Lazy.force

  open B.B

  let pmt_test x1 x2 r1 r2 = 
    let l1 = float x1 /. float r1 
    and l2 = float x2 /. float r2 in 
    let x = 
      if l1 < l2 then Pmt.log_poisson_margin x1 x2 r1 r2
      else Pmt.log_poisson_margin x2 x1 r2 r1
    in -. x
      
  let read_maps = 
    Sle.hfun_make
      (fun x -> lazy (Sam.to_chr_map (Chipseq.bowtie_wodup x))#value)
      Chipseq.samples
      
  let library_size = lazy (
    Sle.hfun_make (fun x -> (Sam.nbmappings (Chipseq.bowtie_wodup x))#value) Chipseq.samples
  )

  let coverage x loc = Chr_map.nbregions_in (! (x --> read_maps)) loc

  let rpkm x loc = 
    float (coverage x loc) /. float (x --> !library_size) *. 1e6 /. (loc |> Location.length |> float) *. 1000.
      
  let pvalue x y loc =
    pmt_test (coverage x loc) (coverage y loc) (x --> !library_size) (y --> !library_size)
end

let all_regions = Target.V.make 
  (object
     method id = "Chipseq_regions.all_regions"
     method deps = [] ++ islands
     method build = 
       let module A = Annotation(struct end) in
       islands#value 
       |> List.map (location_of_island radius)
       |> List.map (fun loc ->
            { loc = loc ;
	      coverage = Sle.hfun_make (fun x -> A.coverage x loc) samples ;
	      rpkm = Sle.hfun_make (fun x -> A.rpkm x loc) samples ;
	      pvalue = Sle.hfun_make (fun (x,y) -> A.pvalue x y loc) paired_samples })
    end)


let bed_of_filter id f = Target.F.make 
  (object
     method id = sprintf "Chipseq_regions.bed_of_filter[%s]" id
     method deps = [] ++ all_regions
     method build path = 
       List.enum all_regions#value 
       // f
       |> Enum.map ((fun x -> x.loc) |- Bed.unparser |- Tsv.string_of_line)
       |> File.write_lines path
     method ty = Tsv.({ has_header = false ; parse = new Bed.base_parser })
   end)


let gerard_filter r = 
  List.exists
    (fun (rar,rxr) -> (rar,`Input_1) --> r.pvalue >= 7. && (rxr,`Input_1) --> r.pvalue >= 7.)
    [ `F9_WT_panRAR_1, `F9_WT_panRXR_1 ;  
      `F9_ATRA_panRAR_1, `F9_ATRA2_panRXR_1 ;
      `F9_ATRA24_panRAR_1, `F9_ATRA24_panRXR_1 ;
      `F9_ATRA48_panRAR_1, `F9_ATRA48_panRXR_1 ]

let selected_regions () = all_regions#value |> List.filter gerard_filter

let selected_regions_bed = 
  bed_of_filter "Chipseq_regions.gerard_filter[r1]" gerard_filter

let gerard_filter_for_time t r = 
  let rar, rxr = match t with 
      `t0  -> `F9_WT_panRAR_1, `F9_WT_panRXR_1
    | `t2  -> `F9_ATRA_panRAR_1, `F9_ATRA2_panRXR_1
    | `t24 -> `F9_ATRA24_panRAR_1, `F9_ATRA24_panRXR_1
    | `t48 -> `F9_ATRA48_panRAR_1, `F9_ATRA48_panRXR_1
  in
  (rar,`Input_1) --> r.pvalue >= 7. && (rxr,`Input_1) --> r.pvalue >= 7.

let string_of_time = function
    `t0  -> "t0"
  | `t2  -> "t2"
  | `t24 -> "t24"
  | `t48 -> "t48"

let selected_regions_bed_for_time t = 
  bed_of_filter 
    (sprintf "Chipseq_regions.gerard_filter_for_time[r1,%s]" (string_of_time t)) 
    (gerard_filter_for_time t)
