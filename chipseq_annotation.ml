open Batteries
open Genome
open Oregon

module Make(P : struct end) = struct
  type annotation = { 
    label : string ;
    f : Location.t -> float
  }

  let ( ! ) = Lazy.force

  let promoter_sel, exon_sel, intron_sel, intragenic_sel = 
    (Ensembl_gff.selection_target B.B.Transcriptome.ensembl_transcripts)#value

  let annotation_of_selection name sel = {
    name ;
    f = RarGenome.Selection.qinclusion sel
  }

  let promoter = annotation_of_selection promoter_sel
  let exon = annotation_of_selection exon_sel
  let intron = annotation_of_selection intron_sel
  let intragenic = annotation_of_selection intragenic_sel

  let tsv_output annotations locz fn = Enum.(locz
    |> map (fun loc ->
	      Array.map (fun x -> x.f loc) annotations)
    |> append 
	(singleton 
	   (Array.map (fun x -> x.name) annotations))
    |> Tsv.unparse_line 
    |> File.write_lines fn
      
end
