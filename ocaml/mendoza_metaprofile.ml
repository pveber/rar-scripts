open Batteries
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

let conditions = [ 
  `Mendoza_metaprofile_RARg;
  `Mendoza_metaprofile_RXRa;
]

let paired_conditions = [ `Mendoza_metaprofile_RARg, `Mendoza_metaprofile_RXRa ]

let all_of_them = 
  Binding_loci.of_macs_targets
    ~radius:250
    (List.map B.Chipseq.macs_peaks conditions)

let gerard_selection =
  let chipseq_annotation =
    Chipseq_annotation.make
      (List.map (fun x -> x, `Mendoza_input) conditions)
      all_of_them
  in V.make
  (object
    method id = "Mendoza_metaprofile.gerard_selection[r1]"
    method deps = [] ++ chipseq_annotation
    method build =
      let ann = chipseq_annotation#value in
      Array.filteri
	(fun i loc ->
	  let pval = ann.(i).Chipseq_annotation.pvalue in
	  List.exists
	    (fun (rar,rxr) -> pval $ rar >= 3. && pval $ rxr >= 3.)
	    paired_conditions)
	all_of_them#value
   end)
  

