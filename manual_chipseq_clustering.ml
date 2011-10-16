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
    (fun (rar,rxr) -> pval x rar >= 7. && pval x rxr >= 7.) 
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
    ~id:"Gerard's selection" ~version:"r4"  
    (fun () -> Array.filter gerard_selection all_regions#value)

type profile = int list * bool list

let comparison_pairs = [
  `F9_WT_panRAR_1, `F9_ATRA_panRAR_1 ; 
  `F9_ATRA_panRAR_1, `F9_ATRA24_panRAR_1 ;
  `F9_ATRA24_panRAR_1, `F9_ATRA48_panRAR_1 ;
  `F9_WT_panRAR_1, `F9_ATRA24_panRAR_1 ;
  `F9_WT_panRAR_1, `F9_ATRA48_panRAR_1 ;
  `F9_ATRA_panRAR_1, `F9_ATRA48_panRAR_1 ;

  `F9_WT_panRXR_1, `F9_ATRA2_panRXR_1 ; 
  `F9_ATRA2_panRXR_1, `F9_ATRA24_panRXR_1 ;
  `F9_ATRA24_panRXR_1, `F9_ATRA48_panRXR_1 ; 
  `F9_WT_panRXR_1, `F9_ATRA24_panRXR_1 ;
  `F9_WT_panRXR_1, `F9_ATRA48_panRXR_1 ;
  `F9_ATRA2_panRXR_1, `F9_ATRA48_panRXR_1 ;
]

let profile_of_region r = 
  List.map (fun (c1, c2) -> compare (intens r c1) (intens r c2)) comparison_pairs,
  List.map (fun c -> pval r c = 0.) Selected_chipseq_regions.conditions

let simple_profile r = 
  List.map (fun (c1, c2) -> compare (intens r c1) (intens r c2)) (List.take 3 comparison_pairs)



let state r x = 
  pval r x > 6.

let state_profile r = [| 
  state r `F9_WT_panRAR_1 ;
  state r `F9_ATRA_panRAR_1 ;
  state r `F9_ATRA24_panRAR_1 ;
  state r `F9_ATRA48_panRAR_1 
|]
  
  
let distribute profile = 
  let regions = Array.filter gerard_selection all_regions#value in
  let accu = Accu.create [] (fun x xs -> x :: xs) in
  Array.iter 
    (fun r -> 
       let p = profile r in 
       Accu.add accu p r)
    regions ;
  List.sort 
    ~cmp:(fun (_,x) (_,y) -> compare (List.length y) (List.length x))
    (List.of_enum (Accu.enum accu))

let cluster_bed f tag n = 
  bed_of_cluster
    ~id:(sprintf "Distribute algorithm[%d]" n) ~version:tag
    (fun () -> snd (List.nth (distribute simple_profile) n) |> Array.of_list)


let control_sequences = 
  Ucsc.control_sequences `mm9 (Ucsc.fasta_of_bed `mm9 gerard_selected)
  |> Fasta.enumerate

let control_enrichment = 
  let ctrlctrlseq = Fasta.enumerate (Fasta.shuffle control_sequences) in
  Motif_library.rank Selected_motifs.value#value control_sequences ctrlctrlseq

let motif_rank bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  Motif_library.rank Selected_motifs.value#value peakseq ctrlseq

let jaspar_motif_rank bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  Motif_library.(rank (of_jaspar_collection Jaspar.core) peakseq ctrlseq)


let meme bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  let prior = MemeSuite.psp_gen ~revcomp:true peakseq ctrlseq in
  MemeSuite.meme ~minw:6 ~maxw:20 ~nmotifs:10 ~revcomp:true ~mode:`zoops ~psp:prior peakseq

  
