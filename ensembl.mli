open Guizmin

(** Transcript annotation downloaded from ftp://ftp.ensembl.org/pub/current_gtf/ 
    with a minor change: chromosome names are changed to chr<n> instead of <n> *)
val transcripts : [`mouse] -> gtf file
