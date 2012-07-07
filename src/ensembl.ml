open Batteries
open Printf
open Biocaml
open Guizmin_bioinfo.MBSchema

type species = [
| `mus_musculus
]


(* acronym of the lab where the species was sequenced *)
let lab_label_of_species = function
| `mus_musculus -> "NCBIM37"

let string_of_species = function
| `mus_musculus -> "mus_musculus"

let fetch_gtf ~release ~species ~path =
  let url = 
    sprintf "ftp://ftp.ensembl.org/pub/release-%d/gtf/%s/%s.%s.%d.gtf.gz"
      release (string_of_species species) 
      (String.capitalize (string_of_species species))
      (lab_label_of_species species) release
  in
  Sle.with_tmp_filename (fun tmp ->
    Sle.sh "wget -O %s %s" tmp url ;
    Sle.sh "zcat %s > %s" tmp path)

let stranded_location_of_gtf_row row = Gtf.(
  row.loc,
  match row.strand with
    | Some `plus -> `Sense
    | Some `minus -> `Antisense
    | _ -> raise (Invalid_argument "Ensembl_gff.stranded_location_of_row")
)

let attr x row =
  List.assoc x row.Gtf.attr

let exon_number = attr "exon_number"
let gene_id = attr "gene_id"
let transcript_id = attr "transcript_id"

let transcripts_enum_aux rows = 
  match List.of_enum rows with
  | [] -> raise (Invalid_argument "Ensembl.transcripts_enum")
  | h :: _ as l ->
      { Transcript.id = transcript_id h ;
        gene_id = gene_id h ;
        strand = (match h.Gtf.strand with
                  | Some `plus -> `Sense
                  | Some `minus -> `Antisense
                  | _ -> assert false) ;
        exons = (
          let compare_exons r1 r2 = 
            compare (exon_number r1) (exon_number r2) in
          let sorted_rows = List.sort compare_exons l in
          List.map (fun r -> r.Gtf.loc) sorted_rows
        ) }


let transcripts_enum gtf = Gtf.(
  (gtf // (fun row -> row.kind = `exon))
  |> Enum.group (fun row -> List.assoc "transcript_id" row.attr)
  |> Enum.map transcripts_enum_aux
)

let tss_map_of_gtf gtf = Gtf.(
  transcripts_enum gtf
  /@ (fun x -> Transcript.tss x, x)
  |> GenomeMap.LMap.of_enum
)


















