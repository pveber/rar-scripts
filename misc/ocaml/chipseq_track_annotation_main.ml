open Batteries
open Biocaml

let bedfn = Sys.argv.(1)
let output = Sys.argv.(2)

let _ =
  Chipseq_track_annotation.(
    tsv
      (all_tracks ())
      (Oregon.Tsv.enum'
         ~has_header:true
         (fun l -> Genome.Location.make l.(0) (int_of_string l.(1)) (int_of_string l.(2)))
         bedfn)
      output
  )










