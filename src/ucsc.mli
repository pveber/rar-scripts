open Batteries

type genome = [ `mm9 | `hg18 ]
val string_of_genome : genome -> string

val fetch_chrom_info : genome -> string -> unit
val chrom_info_of_file : string -> (string * int) Enum.t










