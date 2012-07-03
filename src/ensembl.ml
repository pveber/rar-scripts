open Batteries
open Printf
open Biocaml

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

let tss_map_of_gtf gtf = Gtf.(
  gtf
  // (fun row -> row.kind = `exon && List.assoc "exon_number" row.attr = "1")
  /@ (fun row -> row.loc, (List.assoc "gene_id" row.attr, 
                           List.assoc "transcript_id" row.attr,
                           match row.strand with
                           | Some `plus -> `Sense
                           | Some `minus -> `Antisense
                           | _ -> assert false
  ))
  /@ (fun (loc, ((_,_,dir) as transcript)) -> 
    let loc = Location.upstream ~up:1 ~down:0 dir loc in 
    loc, transcript)
  |> GenomeMap.LMap.of_enum
)


















