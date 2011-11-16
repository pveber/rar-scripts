let promoter_selection = 
  Genome.Selection.of_locations (Ensembl_gff.promoters ~up:1000 Plan.ensembl_gff)
