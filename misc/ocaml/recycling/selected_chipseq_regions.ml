open Batteries
open Oregon
open Genome
open Sle.Infix
open Printf
open Target.Infix

(*
stupid

> sum((d$F9.WT.panRAR.1.input.1..pval. <= -7. & d$F9.WT.panRXR.1.input.1..pval. <= -7.) | (d$F9.ATRA.panRAR.1.input.1..pval. <= -7. & d$F9.ATRA2.panRXR.1.input.1..pval. <= -7) | (d$F9.ATRA24.panRAR.1.input.1..pval. <= -7 & d$F9.ATRA24.panRXR.1.input.1..pval. <= -7) | (d$F9.ATRA48.panRAR.1.input.1..pval. <= -7. & d$F9.ATRA48.panRX)
d$F9.ATRA48.panRXR.1.input.1..pval.  d$F9.ATRA48.panRXR.1...reads.        d$F9.ATRA48.panRXR.1..RPKM.
> sum((d$F9.WT.panRAR.1.input.1..pval. <= -7. & d$F9.WT.panRXR.1.input.1..pval. <= -7.) | (d$F9.ATRA.panRAR.1.input.1..pval. <= -7. & d$F9.ATRA2.panRXR.1.input.1..pval. <= -7) | (d$F9.ATRA24.panRAR.1.input.1..pval. <= -7 & d$F9.ATRA24.panRXR.1.input.1..pval. <= -7) | (d$F9.ATRA48.panRAR.1.input.1..pval. <= -7. & d$F9.ATRA48.panRXR.1.input.1..pval. <= -7))

*)

let conditions = [ 
  `F9_WT_panRAR_1 ; 
  `F9_ATRA_panRAR_1 ;
  `F9_ATRA24_panRAR_1 ;
  `F9_ATRA48_panRAR_1 ;
  `F9_WT_panRXR_1 ; 
  `F9_ATRA2_panRXR_1 ;
  `F9_ATRA24_panRXR_1 ; 
  `F9_ATRA48_panRXR_1
]

let paired_conditions = [
  `F9_WT_panRAR_1, `F9_WT_panRXR_1 ;  
  `F9_ATRA_panRAR_1, `F9_ATRA2_panRXR_1 ;
  `F9_ATRA24_panRAR_1, `F9_ATRA24_panRXR_1 ;
  `F9_ATRA48_panRAR_1, `F9_ATRA48_panRXR_1
]

module Region = struct 
  type t = {
    loc : Location.t ;
    pval : float array ;
    intens : float array
  }

  let index = function
      `F9_WT_panRAR_1 -> 0
    | `F9_ATRA_panRAR_1 -> 1
    | `F9_ATRA24_panRAR_1 -> 2
    | `F9_ATRA48_panRAR_1 -> 3
    | `F9_WT_panRXR_1 -> 4
    | `F9_ATRA2_panRXR_1 -> 5
    | `F9_ATRA24_panRXR_1 -> 6
    | `F9_ATRA48_panRXR_1 -> 7

  let pval r x = r.pval.(index x)
  let intens r x = r.intens.(index x)

end
    
module L = Labs.Make(B.B)
open B.B

let regions = L.binding_loci conditions
let regions_bed = L.binding_loci_tsv conditions

  
let pmt_control ~treatment ~control ~treatment_size ~control_size = 
  Pmt.log_poisson_margin control treatment control_size treatment_size
    

open Region


let binding_loci_pmt n = 
  let read_counts =
    Sle.hfun_make
      (fun x -> (Ngs_utils.read_counts_of_location_file (Chipseq.bowtie_wodup x) regions_bed)#value)
      Chipseq.samples

  and library_size = Sle.hfun_make (fun x -> (Sam.nbmappings (Chipseq.bowtie_wodup x))#value) Chipseq.samples

  and peak_islands = regions#value |> Array.of_list
  in
  let regions = Array.map (fun l -> let x,_,_ = L.peak_average l in Location.(make x.chr (x.st - n) (x.ed + n))) peak_islands
  in
  let profiles_against_input = 
    Array.mapi
      (fun i _ -> 
	 let control = (`Input_1 --> read_counts).(i)
	 and control_size = `Input_1 --> library_size in
	 Sle.hfun_make 
	   (fun x -> -. pmt_control ~treatment:(x --> read_counts).(i) ~control ~treatment_size:(x --> library_size) ~control_size)
	   conditions)
      regions
  in
  let intensity = 
    Array.mapi
      (fun i _ -> 
	 Sle.hfun_make 
	   (fun x -> float (x --> read_counts).(i) /. float (x --> library_size) *. 1e6)
	   conditions)
      regions 
  in
  let samples = Array.of_list conditions in
  (0 --^ Array.length regions) 
  |> Enum.map (fun i -> Location.(let l = regions.(i) in
				  { loc = l ;
				    pval = Array.map (fun x -> x --> profiles_against_input.(i)) samples ;
				    intens = Array.map (fun x -> x --> intensity.(i)) samples }))
  |> Array.of_enum

let target n = Target.V.make
  (object
     method id = sprintf "Selected_chipseq_regions.target[r1,n=%d]" n
     method deps = [] ++ regions_bed ++ regions
     method build = binding_loci_pmt n
   end)

let _ = Pmt.prime_cache()

let annotate_with_chipseq_scores bed path = 
  let validated_chipseq_samples = 
    [ `F9_ATRA_panRAR_1 ;
      `F9_WT_panRAR_1 ; `F9_WT_panRXR_1 ; `F9_ATRA24_panRAR_1 ;
      `F9_ATRA24_panRXR_1 ; `F9_ATRA48_panRAR_1 ; `F9_ATRA48_panRXR_1 ;
      `F9_ATRA2_panRXR_1 ] in
  let read_counts =
    Sle.hfun_make
      (fun x -> (Ngs_utils.read_counts_of_location_file (Chipseq.bowtie_wodup x) regions_bed)#value)
      Chipseq.samples in
  let library_size = Sle.hfun_make (fun x -> (Sam.nbmappings (Chipseq.bowtie_wodup x))#value) Chipseq.samples 
  in
  let pmt_control ~treatment ~control ~treatment_size ~control_size = 
    Pmt.log_poisson_margin control treatment control_size treatment_size 
  in
  let regions = Tsv.enum bed |> Enum.map (fun x -> x # loc) |> Array.of_enum 
  in
  let profiles_against_input = 
    Array.mapi
      (fun i _ -> 
	 let control = (`Input_1 --> read_counts).(i)
	 and control_size = `Input_1 --> library_size in
	 Sle.hfun_make 
	   (fun x -> -. pmt_control ~treatment:(x --> read_counts).(i) ~control ~treatment_size:(x --> library_size) ~control_size)
	   validated_chipseq_samples)
      regions
  in
  let intensity = 
    Array.mapi
      (fun i _ -> 
	 Sle.hfun_make 
	   (fun x -> float (x --> read_counts).(i) /. float (x --> library_size) *. 1e6)
	   validated_chipseq_samples)
      regions
  in
  let samples = Array.of_list validated_chipseq_samples in
  (0 --^ Array.length regions) 
  |> Enum.map (fun i -> Location.(let l = regions.(i) in
				  Array.concat [
				    [| l.chr ; string_of_int l.st ; string_of_int l.ed |] ;
				    Array.map (fun x -> sprintf "%g" (x --> profiles_against_input.(i))) samples ;
				    Array.map (fun x -> sprintf "%g" (x --> intensity.(i))) samples ]))
  |> Enum.append (Enum.singleton (
		    Array.concat [
		      [| "chr" ; "start" ; "end" |] ;
		      Array.map (fun x -> (Chipseq.string_of_sample x) ^ " pval") samples ;
		      Array.map (fun x -> (Chipseq.string_of_sample x) ^ " intens") samples ] ))
  |> Enum.map Tsv.string_of_line
  |> File.write_lines path
