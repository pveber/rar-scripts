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

let panRAR_region_summits_map () = 
  panRAR_regions ()
  |> Enum.mapi (fun i loc -> Location.center loc, i)
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
      |> Array.of_enum)


(** {5 Closest gene from region} *)
let closest_gene_from_region =
  Backup.file
    "rush/chipseq/regions/closest_gene.tsv"
    (fun path -> 
      (panRAR_regions ()
          |> Region_assoc.closest identity (Backup.get tss_map))
      /@ (fun ((chr_region,r_region) as region, ((_,r_tss),(gene_id,_,dir))) -> 
        let (_,r_region_center) = Location.center region in 
        [| 
          chr_region ; 
          string_of_int r_region.Biocaml_range.lo ;
          string_of_int r_region.Biocaml_range.hi ;
          gene_id ;
          string_of_int (Region_assoc.stranded_range_pos ~from:(r_tss, dir) r_region_center) 
        |])
      /@ Tsv.string_of_row
      |> File.write_lines path)

let gene_region_assoc_score = 
  Backup.file
    "rush/chipseq/regions/gene_region_assoc_score.tsv"
    (fun path ->
      Region_assoc.score (chrom_size ())
        fst (LMap.enum (Backup.get tss_map))
        identity (panRAR_regions ())
      /@ (fun (((_,tss_r),(transcript_id,gene_id,dir)),(region_chr, region_r),score) -> [|
        region_chr ; 
        string_of_int region_r.Biocaml_range.lo ;
        string_of_int region_r.Biocaml_range.hi ;
        transcript_id ; gene_id ;
        string_of_int (Region_assoc.range_pos region_r tss_r) ;
        sprintf "%.2f" score |])
      /@ Tsv.string_of_row
      |> File.write_lines path)


(** {5 Genes, their transcripts and expression} *)

(** Amandine's gene table *)
let rnaseq_table () = 
  Rnaseq_table.of_file "resources/rnaseq/synthese_amandine.tsv"

(** A hashtable that gives for each gene_id the corresponding transcripts and their TSS *)
let gene2transcripts () = 
  LMap.enum (Backup.get tss_map)
  /@ (fun (tss_loc, (gene_id, transcript_id, dir)) -> gene_id, (tss_loc, transcript_id, dir))
  |> Biocaml_accu.relation


(** {5 Regions that fall near one of a gene's TSSs} *)

(* FIXME la table d'Amandine est faite avec une version ultérieure à
   la 63, d'où le try with *)
let genes_with_neighbouring_regions () = 
  let regions_map = panRAR_region_summits_map () in 
  let g2t = Hashtbl.of_enum (gene2transcripts ()) in
  rnaseq_table ()
  //@ (fun row ->
    try
      Some (row,
            (Hashtbl.find g2t row.Rnaseq_table.ensembl_gene_id
                |> List.enum)
            /@ (fun (tss_loc, transcript_id,_) -> 
              Region_assoc.neighbours identity 10000 regions_map tss_loc
                   |> Array.enum)
            |> Enum.concat
            |> Set.of_enum)
    with Not_found -> None)




(** {5 Genes: dependence between modulated and region in 10kb radius of tss}*)
let count p l = 
  List.fold_left 
    (fun accu x -> if p x then succ accu else accu)
    0 l 

let modulated g = Rnaseq_table.(
  (abs_float g.log2_fold_6H >=1. && g.padj_6H <= 0.05) ||
  (abs_float g.log2_fold_12H >=1. && g.padj_12H <= 0.05) ||
  (abs_float g.log2_fold_24H >=1. && g.padj_24H <= 0.05) ||
  (abs_float g.log2_fold_36H >=1. && g.padj_36H <= 0.05) ||
  (abs_float g.log2_fold_48H >=1. && g.padj_48H <= 0.05)
)

let region_has_motif motif_annotation motif i =
  List.exists
    (fun (mot,_,_,_) -> motif = mot)
    motif_annotation.(i)

let gene_has_motif motif_annotation motif (_, regions) =
  Set.exists
    (snd |- region_has_motif motif_annotation motif)
    regions


let test ?(filter = fun _ -> true) p1 p2 xs =
  let open Biocaml_accu in
  let counter = Counter.of_enum (xs // filter /@ (fun x -> p1 x, p2 x)) in
  let pp = get counter (true,  true)
  and pn = get counter (true,  false)
  and np = get counter (false, true)
  and nn = get counter (false, false)
  in 
  printf "(%d %d) (%d %d)\n" pp pn np nn


let gene_test_modulated_vs_region_in_promoter () = 
  test
    
    (fst |- modulated)
    (fun (_,regions) -> Set.cardinal regions > 0)
    (genes_with_neighbouring_regions ())

let gene_test_modulated_vs_motif motif =
  let motif_annotation = Backup.get panRAR_regions_motif_annotation in 
  test
    (fst |- modulated)
    (gene_has_motif motif_annotation motif)
    (genes_with_neighbouring_regions ())


(* let () = *)
(*   gene_test_modulated_vs_region_in_promoter () *)

let () =
  gene_test_modulated_vs_motif (`direct, 0)



let () =
  ignore (Backup.get closest_gene_from_region)



















