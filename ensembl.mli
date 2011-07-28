open Batteries
open Guizmin

module Guizmin_plugin : sig 

  val transcripts : [`mouse] -> gtf file
(** Transcript annotation downloaded from ftp://ftp.ensembl.org/pub/current_gtf/ 
    with a minor change: chromosome names are changed to chr<n> instead of <n> *)

end

val promoters : gtf file -> Genome.location Enum.t
