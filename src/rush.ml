open Batteries

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

let _ = 
  Region_assoc.closest identity (Backup.get tss_map) (List.enum panRAR_peaks)
  |> List.of_enum











