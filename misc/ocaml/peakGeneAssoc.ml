open Batteries
open Biocaml
open Genome
open HFun.Infix

module LMap = Biocaml.GenomeMap.LMap

let closest_ensembl_gene ls = 
  let tssmap = Ensembl_gff.tss B.Transcriptome.ensembl_transcripts in 
  let rank = Ensembl_gff.rank_function B.Transcriptome.ensembl_transcripts in
  ls
  /@ (fun loc -> LMap.closest (S.newloc loc) tssmap)
  /@ (fun ((chr,r),x,d) -> S.oldloc (chr,r), x, d, try rank $ x with _ -> failwith x)

let panRAR_regions_to_gene_tsv path = 
  let regions  = PanRAR_regions.gerard_selection#value in 
  let closest_genes = closest_ensembl_gene (Array.enum regions) |> Array.of_enum in 
  let pos = Array.map2 (fun r (loc,gene_id,_,rank) -> loc, Location.position ~from:r loc, gene_id,rank) regions closest_genes in
  Array.enum pos
  /@ Location.(fun (loc, pos, gene_id,rank) -> [| loc.chr ; string_of_int loc.st ; string_of_int loc.ed ; string_of_int pos ; gene_id ; string_of_int rank|])
  /@ Oregon.Tsv.string_of_line
  |> Enum.append (Enum.singleton "chr\tstart\tend\tposition\tgene_id\tgene_rank")
  |> File.write_lines path
  
  




















