let all_regions_fa = 
  let path = "manual/chipseq/regions/all_regions.bed" in
  (all_regions ()) /@ Bed.unparser /@ Tsv.string_of_line |> File.write_lines path ;	(* UGLY ! *)		
  let bed = Target.F.input path Tsv.({ has_header = false ; parse = new Bed.base_parser }) in 
  Ucsc.fasta_of_bed `mm9 bed

let format_res l = 
  List.sort ~cmp:(fun (_,_,x,_) (_,_,y,_) -> Pervasives.compare y x) l
  |> List.take 5
  |> List.map (fun (st,ed,x,pv) -> [| string_of_int st ; string_of_int ed ; string_of_int x ; string_of_float pv |])
  |> Array.concat 
  |> Tsv.string_of_line

let tsv_output item fa fn =
  let r = (B.B.TFBS.wapam_fasta_search 1e-2 item fa)#value in
  Array.enum r
  /@ format_res
  |> File.write_lines fn

let make locz prefix = 
  let fa = 

  Array.iter 
    (fun item -> tsv_output item locz (sprintf "%s-%s.bed" prefix (Motif_library.string_of_item item)))
    Selected_motifs.drerir
