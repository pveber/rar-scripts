open Batteries
open Biocaml
(* open Guizmin *)

(*
type item = Exon of string * int
type annotation = item Biocaml.Genome.Annotation.t

let annotation = assert false
*)

let stranded_location_of_row row = Gff.(
  RarGenome.Location.make row.chr (fst row.pos) (snd row.pos),
  match row.strand with
    | Sense -> `Sense
    | Antisense -> `Antisense
    | _ -> raise (Invalid_argument "Ensembl.stranded_location_of_row")
)

let promoters ?(up = 1000) ?(down = 0) path = Gff.(
  List.enum (to_list (of_file ~version:2 path))
  |> Enum.filter 
      (fun row -> try row.feature = "exon" 
		  && get_attribute row "exon_number" <> "1" with _ -> failwith (row_to_string row))
  |> Enum.map stranded_location_of_row
  |> Enum.map (fun (loc,strand) -> RarGenome.Location.upstream ~up ~down strand loc)
)

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
