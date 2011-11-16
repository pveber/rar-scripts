open Batteries

let chipseq_design = 
  List.concat [
    List.map (fun x -> x, `Input_1) PanRAR_regions.conditions ;
    List.map (fun x -> x, `Mendoza_input) Rarg_regions.conditions
  ]
    

let original_panRAR_regions = Oregon.(
  Target.F.input 
    "results/chipseq/regions/panRAR_regions.bed" 
    (Bed.basic_ty (assert false))
)



(** CHIPSEQ *)
let _ = Chipseq_annotation.(
  tsv
    "results/chipseq/annotations/RARg_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       Rarg_regions.gerard_selection)#value
);;

let _ = Chipseq_annotation.(
  tsv
    "results/chipseq/annotations/PanRAR_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       PanRAR_regions.gerard_selection)#value
);;

let _ = Chipseq_annotation.(
  tsv
    "results/chipseq/annotations/original_PanRAR_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       (S.bed_to_value original_panRAR_regions))#value
)

(** CONSERVATION *)
let _ = Conservation_annotation.make 
  (Array.enum (S.bed_to_value original_panRAR_regions)#value)
  "results/chipseq/annotations/original_PanRAR_regions_conservation_annotation.tsv"

let _ = Conservation_annotation.make 
  (Array.enum Rarg_regions.gerard_selection#value)
  "results/chipseq/annotations/RARg_regions_conservation_annotation.tsv"
