open Batteries
(* open Guizmin *)
open Genome
open Oregon
open Target.Infix

(*
type item = Exon of string * int
type annotation = item Biocaml.Genome.Annotation.t

let annotation = assert false
*)

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
