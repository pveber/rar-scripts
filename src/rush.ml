open Batteries
open Printf
module LMap = Biocaml_genomeMap.LMap

let chrom_info_file = Backup.file
  "rush/ucsc/mm9/chrom_info"
  (Ucsc.fetch_chrom_info `mm9)

let chrom_info () = 
  Ucsc.chrom_info_of_file (Backup.get chrom_info_file)
  |> List.of_enum

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

let panRAR_peaks = 
  Bed.enum ~header:true "resources/chipseq/regions/PanRAR_regions.bed"
  |> List.of_enum

let closest_gene_from_peak =
  Backup.file
    "rush/chipseq/peaks/closest_gene.tsv"
    (fun path -> 
      (List.enum panRAR_peaks
          |> Region_assoc.closest identity (Backup.get tss_map))
      /@ (fun ((chr_peak,r_peak), ((_,r_tss),(gene_id,_))) -> [| 
        chr_peak ; 
        string_of_int r_peak.Biocaml_range.lo ;
        string_of_int r_peak.Biocaml_range.hi ;
        gene_id ;
        string_of_int (Region_assoc.range_pos r_peak r_tss) |])
      /@ Tsv.string_of_row
      |> File.write_lines path)

let _ = 
  Region_assoc.score (chrom_info ())
    fst ((LMap.enum (Backup.get tss_map)) // (fun ((chr,_),_) -> chr = "chr10"))
    identity (List.enum panRAR_peaks)
  |> Enum.iter 
      (fun (_,_,score) -> printf "%g\n%!" score)
















