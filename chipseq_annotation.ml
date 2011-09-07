open Batteries
open Printf
open Genome
open Oregon
open Sle.Infix

module Make(P : sig end) = struct
  type annotation = { 
    label : string ;
    f : Location.t -> string
  }

  let ( ! ) = Lazy.force

  let location = {
    label = "Location" ;
    f = Location.to_string
  }


  (*************************** sequence kinds ************************************)
  let promoter_sel, exon_sel, intron_sel, intragenic_sel = 
    (Ensembl_gff.selection_target B.B.Transcriptome.ensembl_transcripts)#value

  let annotation_of_selection label sel = {
    label ;
    f = (fun loc -> sprintf "%.2f" (RarGenome.Selection.qinclusion sel loc))
  }



  let promoter = annotation_of_selection "Promoter region (1kb)" promoter_sel
  let exon = annotation_of_selection "Exonic region" exon_sel
  let intron = annotation_of_selection "Intronic region" intron_sel
  let intragenic = annotation_of_selection "Intragenic region" intragenic_sel




  (*************************** chipseq coverage ************************************)
  module Chipseq_coverage_internals = struct
    open B.B
    open Sle.Infix

    let _ = Pmt.prime_cache()

    let pmt_test x1 x2 r1 r2 = 
      let l1 = float x1 /. float r1 
      and l2 = float x2 /. float r2 in 
      if l1 < l2 then Pmt.log_poisson_margin x1 x2 r1 r2
      else Pmt.log_poisson_margin x2 x1 r2 r1

    let read_maps = 
      Sle.hfun_make
	(fun x -> lazy (Sam.to_chr_map (Chipseq.bowtie_wodup x))#value)
	Chipseq.samples

    let library_size = lazy (
      Sle.hfun_make (fun x -> (Sam.nbmappings (Chipseq.bowtie_wodup x))#value) Chipseq.samples
    )

    let coverage x loc = Chr_map.nbregions_in (! (x --> read_maps)) loc

    let rpkm x loc = 
      float (coverage x loc) /. float (x --> !library_size) *. 1e6 /. (loc |> Location.length |> float) *. 1000.

    let pvalue x y loc =
      pmt_test (coverage x loc) (coverage y loc) (x --> !library_size) (y --> !library_size)
  end

  let coverage s = { 
    label = sprintf "%s (#reads)" (B.B.Chipseq.string_of_sample s) ; 
    f = fun loc -> Chipseq_coverage_internals.coverage s loc |> string_of_int
  }

  let rpkm s = { 
    label = sprintf "%s (RPKM)" (B.B.Chipseq.string_of_sample s) ; 
    f = fun loc -> Chipseq_coverage_internals.rpkm s loc |> string_of_float
  }

  let pvalue x y = { 
    label = sprintf "%s/%s (pval)" (B.B.Chipseq.string_of_sample x) (B.B.Chipseq.string_of_sample y) ; 
    f = fun loc -> Chipseq_coverage_internals.pvalue x y loc |> string_of_float
  }

  let validated_chipseq_samples = [| 
    `F9_WT_panRAR_1 ; `F9_WT_panRXR_1 ; 
    `F9_ATRA_panRAR_1 ; `F9_ATRA2_panRXR_1 ;
    `F9_ATRA24_panRAR_1 ; `F9_ATRA24_panRXR_1 ; 
    `F9_ATRA48_panRAR_1 ; `F9_ATRA48_panRXR_1 ;
    `Input_1 |]

  let comparison_pairs = [| 
    `F9_WT_panRAR_1, `Input_1 ; 
    `F9_WT_panRXR_1, `Input_1 ; 
    `F9_ATRA_panRAR_1, `Input_1 ; 
    `F9_ATRA2_panRXR_1, `Input_1 ;
    `F9_ATRA24_panRAR_1, `Input_1 ; 
    `F9_ATRA24_panRXR_1, `Input_1 ;
    `F9_ATRA48_panRAR_1, `Input_1 ; 
    `F9_ATRA48_panRXR_1, `Input_1 ;
    `F9_WT_panRAR_1, `F9_ATRA_panRAR_1 ;
    `F9_WT_panRXR_1, `F9_ATRA2_panRXR_1 ;
    `F9_ATRA_panRAR_1, `F9_ATRA24_panRAR_1 ;
    `F9_ATRA2_panRXR_1, `F9_ATRA24_panRXR_1 ;
    `F9_ATRA24_panRAR_1, `F9_ATRA48_panRAR_1 ;
    `F9_ATRA24_panRXR_1, `F9_ATRA48_panRXR_1 ;
    `F9_WT_panRAR_1, `F9_ATRA24_panRAR_1 ;
    `F9_WT_panRXR_1, `F9_ATRA24_panRXR_1 ;
    `F9_WT_panRAR_1, `F9_ATRA48_panRAR_1 ;
    `F9_WT_panRXR_1, `F9_ATRA48_panRXR_1 ;
    `F9_ATRA_panRAR_1, `F9_ATRA48_panRAR_1 ;
    `F9_ATRA2_panRXR_1, `F9_ATRA48_panRXR_1 ;
  |]










  (*************************** neighboring genes ************************************)
  module Neighboring_genes_internals = struct
    let tss_map = lazy (
      let tss_of_gene = Ucsc.promoter_of_gene `mm9 ~upstream:0 ~downstream:0 in
      let transcripts = Tsv.enum (Ucsc.ensGene `mm9) |> Array.of_enum in
      let t2g = Sle.hrel_of_enum (fun g -> g#transcript, g#gene) (Tsv.enum (Ucsc.ensGtp `mm9)) in
      let tss = Array.map (fun t -> tss_of_gene t, List.hd (t#name --> t2g)) transcripts in 
      Chr_map.of_array tss
    )

    let filter_already_seen f a =
      let accu = ref PSet.empty in
      ((Array.enum a) // (fun x ->
			    let y = f x in 
			    let b = PSet.mem y accu.contents in
			    accu := PSet.add y accu.contents ; not b)) |>
      Array.of_enum
			   


    let k_closest_genes k locs = 
      let r = filter_already_seen (fun (_,g,_) -> g) locs in 
      if Array.length r >= k then Some (Array.sub r 0 k) else None

    let kth_closest_gene k loc =
      (Chr_map.closest_regions (k_closest_genes k) loc !tss_map).(k - 1)
  end

  let kth_closest_tss_gene k = {
    label = sprintf "Gene with closest TSS #%d" k ;
    f = (fun loc -> Neighboring_genes_internals.(
	   try Tuple3.second (kth_closest_gene k loc)
	   with Not_found -> "--"))
  }

  let kth_closest_tss_gene_distance k = {
    label = sprintf "Distance of gene with closest TSS #%d" k ;
    f = (fun loc -> Neighboring_genes_internals.(
	   try Location.position ~from:loc (Tuple3.first (kth_closest_gene k loc)) |> string_of_int
	   with Not_found -> "--"))
  }



  let tsv_output annotations locz fn = Enum.(locz
    |> map (fun loc ->
	      Array.map (fun x -> x.f loc) annotations)
    |> append 
	(singleton 
	   (Array.map (fun x -> x.label) annotations))
    |> map Tsv.string_of_line 
    |> File.write_lines fn
  )


  let all = Array.concat [
    [| location ; promoter ; exon ; intron ; intragenic |] ;
    Array.map coverage validated_chipseq_samples ;
    Array.map rpkm validated_chipseq_samples ;
    Array.map (fun (x,y) -> pvalue x y) comparison_pairs ;
    Array.concat (List.init 4 (fun i -> [| kth_closest_tss_gene (i + 1) ; kth_closest_tss_gene_distance (i + 1) |])) ;
  ]
      
end
