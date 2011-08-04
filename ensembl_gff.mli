open Batteries

(*
open Guizmin

module Guizmin_plugin : sig
  type gff 
    
  val gff : [`mouse] -> gff file
    (** Transcript annotation downloaded from ftp://ftp.ensembl.org/pub/current_gtf/ 
	with a minor change: chromosome names are changed to chr<n> instead of <n> *)
end
*)

(*  
type strand = Biocaml.Gff.strand
type item = Exon of string * int * strand | Other
type annotation = item Genome.Annotation.t
    
val annotation : Guizmin.gff file -> annotation
*)
val promoters : ?up:int -> ?down:int -> string -> RarGenome.Location.t Enum.t

