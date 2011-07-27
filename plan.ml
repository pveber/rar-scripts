open Guizmin

let justin_clustering : unit dir = 
  let archive = Unix_utils.unvac ~host:"titan.ibisc.fr" ~adr:"vac:bec7ce44a0702d261c56f88e31ec2dd6a1128f0d" in
  alias
    (partof archive "clustering")
    [ "justin" ; "clustering" ]

let ensembl_gtf = Ensembl.transcripts `mouse

