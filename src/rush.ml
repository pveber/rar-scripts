let gtf = Result.file 
  "rush/ensembl/gtf/mus_musculus.63.gtf" 
  (fun path -> Ensembl.fetch_gtf ~release:63 ~species:`mus_musculus ~path)

let _ = Result.get gtf




















