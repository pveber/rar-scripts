open Batteries
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

let conditions = [ 
  `Mendoza_ATRA2_RARg;
  `Mendoza_ATRA6_RARg;
  `Mendoza_ATRA24_RARg;
  `Mendoza_ATRA48_RARg;
  `Mendoza_EtOH48_RARg;
  `Mendoza_ATRA2_RXRa;
  `Mendoza_ATRA6_RXRa;
  `Mendoza_ATRA24_RXRa;
  `Mendoza_ATRA48_RXRa;
  `Mendoza_EtOH48_RXRa;
]

let paired_conditions = [
  `Mendoza_ATRA2_RARg,   `Mendoza_ATRA2_RXRa;
  `Mendoza_ATRA6_RARg,   `Mendoza_ATRA6_RXRa;
  `Mendoza_ATRA24_RARg,   `Mendoza_ATRA24_RXRa;
  `Mendoza_ATRA48_RARg,   `Mendoza_ATRA48_RXRa;
  `Mendoza_EtOH48_RARg,   `Mendoza_EtOH48_RXRa;
]

let all_of_them = 
  Binding_loci.of_macs_targets
    ~radius:250
    (List.map B.Chipseq.macs_peaks conditions)


(* let gerard_selection =  *)
(*   let chipseq_annotation =  *)
(*     Chipseq_annotation.make *)
(*       (List.map (fun x -> x, `Input_1) conditions) *)
(*       all_of_them *)
(*   in V.make  *)
(*   (object *)
(*     method id = "RARg_regions.gerard_selection" *)
(*     method deps = [] ++ chipseq_annotation *)
(*     method build =  *)
(*       let ann = chipseq_annotation#value in *)
(*       Array.filteri *)
(* 	(fun i loc ->  *)
(* 	  let pval = ann.(i).Chipseq_annotation.pvalue in *)
(* 	  List.exists  *)
(* 	    (fun (rar,rxr) -> pval $ rar >= 7. && pval $ rxr >= 7.) *)
(* 	    paired_conditions) *)
(* 	all_of_them#value *)
(*    end) *)
  

