open Batteries
open Printf
open Selected_chipseq_regions.Region
open Oregon
open Target.Infix

let conj l x = List.for_all (fun f -> f x) l

let peak_radius = 250
let all_regions = Selected_chipseq_regions.target peak_radius

let gerard_selection x = 
  List.exists 
    (fun (rar,rxr) -> pval x rar > 6. && pval x rxr > 5.) 
    Selected_chipseq_regions.paired_conditions

let bed_of_cluster ~id ~version f = Target.F.make 
  (object
     method id = sprintf "Manual_chipseq_clustering.bed_of_cluster[%s,%s]" id version
     method deps = [] ++ all_regions
     method build path = 
       Array.enum (f ()) 
       |> Enum.map ((fun x -> x.loc) |- Bed.unparser |- Tsv.string_of_line)
       |> File.write_lines path
     method ty = Tsv.({ has_header = false ; parse = new Bed.base_parser })
   end)

let gerard_selected = 
  bed_of_cluster 
    ~id:"Gerard's selection" ~version:"r1"  
    (fun () -> Array.filter gerard_selection all_regions#value)

let control_sequences = 
  Ucsc.control_sequences `mm9 (Ucsc.fasta_of_bed `mm9 gerard_selected)

let motif_rank bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  Motif_library.rank Selected_motifs.value#value peakseq ctrlseq

