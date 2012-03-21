open Batteries
open Biocaml
open Genome

module LMap = Biocaml.GenomeMap.LMap

let closest_ensembl_gene ls = 
  let tssmap = Ensembl_gff.tss B.Transcriptome.ensembl_transcripts in 
  ls
  /@ (fun loc -> LMap.closest (S.newloc loc) tssmap)
  /@ (fun ((chr,r),x,d) -> S.oldloc (chr,r), x, d)

let panRAR_regions_to_gene_tsv path = 
  let regions  = PanRAR_regions.gerard_selection#value in 
  let closest_genes = closest_ensembl_gene (Array.enum regions) |> Array.of_enum in 
  let pos = Array.map2 (fun r (loc,gene_id,_) -> loc, Location.position ~from:r loc, gene_id) regions closest_genes in
  Array.enum pos
  /@ Location.(fun (loc, pos, gene_id) -> [| loc.chr ; string_of_int loc.st ; string_of_int loc.ed ; string_of_int pos ; gene_id |])
  /@ Oregon.Tsv.string_of_line
  |> File.write_lines path
  
  




















