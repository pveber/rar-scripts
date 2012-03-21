open Batteries
open Genome
open Oregon
open Target
open Target.Infix

module LMap = Biocaml_genomeMap.LMap

let stranded_location_of_row row = Gtf.(
  row.loc,
  match row.strand with
    | Some `plus -> `Sense
    | Some `minus -> `Antisense
    | _ -> raise (Invalid_argument "Ensembl_gff.stranded_location_of_row")
)

let promoters ?(up = 1000) ?(down = 0) gff = Gtf.(Tsv.enum gff
  |> Enum.filter 
      (fun row -> row.kind = `exon && List.assoc "exon_number" row.attr = "1")
  |> Enum.map stranded_location_of_row
  |> Enum.map (fun (loc,strand) -> RarGenome.location_upstream ~up ~down strand loc)
  |> RarGenome.Selection.of_locations
)

let tss gff = Gtf.(
  Tsv.enum gff
  // (fun row -> row.kind = `exon && List.assoc "exon_number" row.attr = "1")
  /@ (fun row -> stranded_location_of_row row, List.assoc "gene_id" row.attr)
  /@ (fun ((loc,strand), gene_id) -> RarGenome.location_upstream ~up:1 ~down:0 strand loc, gene_id)
  /@ (fun (loc,gene_id) -> Location.(loc.chr, Biocaml_range.make loc.st loc.ed), gene_id)
  |> LMap.of_enum
)

let exons gff = Gtf.(Tsv.enum gff
  |> Enum.filter_map (fun row -> if row.kind = `exon then Some row.loc else None)
  |> RarGenome.Selection.of_locations
)

let intragenic gff = Gtf.(
  let merge l1 l2 = Location.(
    make l1.chr (min l1.st l2.st) (max l1.ed l2.ed)
  ) in
  Tsv.enum gff 
  |> Enum.filter (fun row -> row.kind = `exon)
  |> Enum.group (fun row -> List.assoc "transcript_id" row.attr)
  |> Enum.map (Enum.map (fun row -> row.loc))
  |> Enum.map (Enum.reduce merge)
  |> RarGenome.Selection.of_locations
)

let introns gff = 
  RarGenome.Selection.diff 
    (intragenic gff)
    (exons gff)

let selection_target gff = Target.V.make
  (object
     method id = "Ensembl_gff.selection_target[r3]"
     method deps = [] ++ gff
     method build =
       (promoters ~up:1000 gff,
	exons gff,
	introns gff,
	intragenic gff)
   end)

let convex_hull_of_transcripts ts =
  Location.convex_hull
    (List.map
       (fun (_,t) -> Location.convex_hull (List.map (fun x -> x.Gtf.loc) t))
       ts)

let gene_from_gff gid transcripts = 
  assert (transcripts <> []) ;
  let _, t1 = List.hd transcripts in
  assert (t1 <> []) ;
  let symbol = List.assoc "gene_name" (List.hd t1).Gtf.attr in 
  Gene.({ 
    id = gid ; symbol ; 
    loc = convex_hull_of_transcripts transcripts ;
    transcripts = 
      List.map
	Transcript2.(fun (_,exons) -> 
		       assert (exons <> []) ;
		       let fst_it = List.hd exons in {
			 id = List.assoc "transcript_id" fst_it.Gtf.attr ;
			 gene_symbol = List.assoc "gene_name" fst_it.Gtf.attr ;
			 description = "" ;
			 exons = List.map (fun it -> it.Gtf.loc) exons ;
			 strand = Option.get fst_it.Gtf.strand
		       })
	transcripts })
    


let genes gff = V.make 
  (object
     method id = "Ensembl_gff.genes[r1]"
     method deps = [] ++ gff
     method build = 
       let transcripts = Gtf.transcripts identity gff in
       (Hashtbl.enum transcripts 
          |> Sle.hrel_of_enum (fun ((gid,tid),transcript) -> gid, (tid, transcript.Transcript.annot))
	  |> Hashtbl.enum) (* transcripts by gid *)
       /@ (fun (gid,ts) -> gene_from_gff gid ts)
       |> Array.of_enum
   end)



(*
module Guizmin_plugin = struct
  type gff 

  let string_of_org = function
  `mouse -> "mouse"
  let transcript_url = function
  `mouse -> "ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.NCBIM37.63.gtf.gz"

  let gff x = target 
    [ "mkdir -p $@_download" ;
      sp "cd $@_download && wget %s" (transcript_url x) ;
      "cd $@_download && gunzip *.gz" ;
      "D=${PWD} && cd $@_download && gawk '{print \"chr\"$0}' *.gtf | sed 's/chrMT/chrM/g' | grep -v 'ENSMUST00000127664' > ${D}/$@" ;
      "rm -rf $@_download" ]
    []
end
*)

