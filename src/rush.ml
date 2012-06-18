open Batteries
open Printf
module LMap = Biocaml_genomeMap.LMap

(** {5 Genes and genome definition} *)

(** Tabular file for chromosome size *)
let chrom_size_file = Backup.file
  "rush/ucsc/mm9/chrom_size"
  (Ucsc.fetch_chrom_size `mm9)

(** [chrom_size ()] returns an array indicating chromosome names and
    length *)
let chrom_size () = 
  Ucsc.chrom_size_of_file (Backup.get chrom_size_file)
  |> Array.of_enum

(** release 63 of ensembl genes *)
let gtf_file = Backup.file 
  "rush/ensembl/63/mus_musculus/gene_set.gtf"
  (fun path -> Ensembl.fetch_gtf ~release:63 ~species:`mus_musculus ~path)

let gtf () = 
  let chr = function
  | "MT" -> "chrM"
  | s -> "chr" ^ s
  in 
  Gtf.enum ~chr (Backup.get gtf_file)

let tss_map = Backup.value
  "rush/ensembl/63/mus_musculus/tss_map.value"
  (fun () -> Ensembl.tss_map_of_gtf (gtf ()))







(** {5 The binding regions} *)

(** As an enumeration *)
let panRAR_regions () = 
  Bed.enum ~header:true "resources/chipseq/regions/PanRAR_regions.bed"

(** As an LMap *)
let panRAR_regions_map () = 
  panRAR_regions ()
  |> Enum.mapi (fun i loc -> loc, i)
  |> LMap.of_enum

(** As a fasta *)
let panRAR_regions_sequences () = 
  (Biocaml_fasta.enum_of_file "resources/chipseq/regions/PanRAR_regions.fa" |> snd)
  /@ snd

(** {5 Motif score annotation for regions} *)
let panRAR_regions_motif_annotation =
  Backup.value
    "rush/chipseq/regions/motif_annotation.value"
    (fun () -> 
      Motif_annotation.of_sequences Nhre.motifs 0.05 (panRAR_regions_sequences ())
      |> List.of_enum)


(** {5 Closest gene from region} *)
let closest_gene_from_region =
  Backup.file
    "rush/chipseq/regions/closest_gene.tsv"
    (fun path -> 
      (panRAR_regions ()
          |> Region_assoc.closest identity (Backup.get tss_map))
      /@ (fun ((chr_region,r_region), ((_,r_tss),(gene_id,_))) -> [| 
        chr_region ; 
        string_of_int r_region.Biocaml_range.lo ;
        string_of_int r_region.Biocaml_range.hi ;
        gene_id ;
        string_of_int (Region_assoc.range_pos r_region r_tss) |])
      /@ Tsv.string_of_row
      |> File.write_lines path)

let gene_region_assoc_score = 
  Backup.file
    "rush/chipseq/regions/gene_region_assoc_score.tsv"
    (fun path ->
      Region_assoc.score (chrom_size ())
        fst ((LMap.enum (Backup.get tss_map)) // (fun ((chr,_),_) -> chr = "chr10"))
        identity (panRAR_regions ())
      /@ (fun (((_,tss_r),(transcript_id,gene_id)),(region_chr, region_r),score) -> [|
        region_chr ; 
        string_of_int region_r.Biocaml_range.lo ;
        string_of_int region_r.Biocaml_range.hi ;
        transcript_id ; gene_id ;
        string_of_int (Region_assoc.range_pos region_r tss_r) ;
        sprintf "%.2f" score |])
      /@ Tsv.string_of_row
      |> File.write_lines path)


(** {5 Closest region from gene} *)

(** Amandine's gene table *)
let rnaseq_table () = 
  Rnaseq_table.of_file "resources/rnaseq/synthese_amandine.tsv"

let gene_neighbouring_regions () = 
  Region_assoc.neighbours
    Rnaseq_table.tss_loc_of_row
    10000
    (panRAR_regions_map ())
    (rnaseq_table ())




(** {5 Genes: dependence between modulated and region in 10kb radius of tss}*)
let count p l = 
  List.fold_left 
    (fun accu x -> if p x then succ accu else accu)
    0 l 

let modulated g = Rnaseq_table.(
  (abs_float g.log2_fold_6H > 2. && g.padj_6H < 0.05) ||
  (abs_float g.log2_fold_12H > 2. && g.padj_12H < 0.05) ||
  (abs_float g.log2_fold_24H > 2. && g.padj_24H < 0.05) ||
  (abs_float g.log2_fold_36H > 2. && g.padj_36H < 0.05) ||
  (abs_float g.log2_fold_48H > 2. && g.padj_48H < 0.05)
)

let gene_test_modulated_vs_region_in_promoter () = 
  let gene_and_neighbours = List.of_enum (gene_neighbouring_regions ()) in 
  let n = List.length gene_and_neighbours 
  and k_mod = 
    count
      (fun (gene, regions) -> modulated gene && regions <> [||])
      gene_and_neighbours
  and k_non_mod = 
    count
      (fun (gene, regions) -> not (modulated gene) && regions <> [||])
      gene_and_neighbours
  and n_mod = 
    count (fun (gene, _) -> modulated gene) gene_and_neighbours
  and n_non_mod = 
    count (fun (gene, _) -> not (modulated gene)) gene_and_neighbours
  in 
  printf "(%d %d) (%d %d) %f %f\n" k_mod n_mod k_non_mod n_non_mod (float k_mod /. float k_non_mod) (float k_non_mod /. float n_non_mod)


let () =
  let gene_and_neighbours = List.of_enum (gene_neighbouring_regions ()) in  
  let open Rnaseq_table in
      List.iter (fun (gene, regions) ->
        printf "%s %s %d %d %s\n%!" gene.gene_symbol gene.chr gene.chr_start gene.chr_end (Location.to_string (tss_loc_of_row gene)) ;
        Array.iter
          (fun (loc,_) -> print_endline (Location.to_string loc))
          regions
      )
        gene_and_neighbours

(*
let () = 
  gene_test_modulated_vs_region_in_promoter ()
*)


(* let () =  *)
(*   ignore (Backup.get closest_gene_from_region) *)



















