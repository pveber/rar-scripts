let chipseq_design = 
  List.concat [
    List.map (fun x -> x, `Input_1) PanRAR_regions.conditions ;
    List.map (fun x -> x, `Mendoza_input) Rarg_regions.conditions
  ]
    
let _ = Chipseq_annotation.(
  tsv
    "RARg_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       Rarg_regions.gerard_selection)#value
);;

let _ = Chipseq_annotation.(
  tsv
    "PanRAR_regions_chipseq_annotation.tsv"
    (make 
       chipseq_design
       PanRAR_regions.gerard_selection)#value
);;
