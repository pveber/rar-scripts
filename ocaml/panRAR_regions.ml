open Batteries
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

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

let all_of_them = 
  Binding_loci.of_macs_targets
    ~radius:250
    (List.map B.Chipseq.macs_peaks conditions)

let gerard_selection = 
  let chipseq_annotation = 
    Chipseq_annotation.make
      (List.map (fun x -> x, `Input_1) conditions)
      all_of_them
  in V.make 
  (object
    method id = "PanRAR_regions.gerard_selection"
    method deps = [] ++ chipseq_annotation
    method build = 
      let (values,_,_) as ann = chipseq_annotation#value in
      Array.filteri
	(fun i loc -> 
	  let pval = values.(i).Chipseq_annotation.pvalue
          and i2 = Chipseq_annotation.i2 ann in
	  List.exists 
	    (fun (rar,rxr) -> pval.(i2 rar `Input_1) >= 7. && pval.(i2 rxr `Input_1) >= 7.)
	    paired_conditions)
	all_of_them#value
   end)
  
















