open Batteries

let chipseq_design = 
  List.concat [
    List.map (fun x -> x, `Input_1) PanRAR_regions.conditions ;
    List.map (fun x -> x, `Mendoza_input) Rarg_regions.conditions ;
    [ `Mendoza_metaprofile_RARg, `Mendoza_input; `Mendoza_metaprofile_RXRa, `Mendoza_input ]
  ]
    

let original_panRAR_regions = Oregon.(
  Target.F.input 
    "resources/chipseq/regions/original_PanRAR_regions.bed" 
    (Bed.basic_ty false)
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

let _ = Chipseq_annotation.(
  tsv
    "results/chipseq/annotations/metaprofile_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       Mendoza_metaprofile.gerard_selection)#value
);;

(** CONSERVATION *)
let _ = Conservation_annotation.make 
  (Array.enum (S.bed_to_value original_panRAR_regions)#value)
  "results/chipseq/annotations/original_PanRAR_regions_conservation_annotation.tsv"

let _ = Conservation_annotation.make 
  (Array.enum Rarg_regions.gerard_selection#value)
  "results/chipseq/annotations/RARg_regions_conservation_annotation.tsv"


(** RNASEQ **)
let () = Rnaseq_annotation.(
  tsv 
    "results/rnaseq/ensembl_rnaseq_annotation.tsv" 
    (make 
       B.Setting.RnaseqData.design 
       (Ensembl_gff.genes B.Transcriptome.ensembl_transcripts))
)
