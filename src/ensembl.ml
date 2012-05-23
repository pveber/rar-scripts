open Batteries
open Printf

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






















