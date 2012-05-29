open Batteries

type genome = [ `mm9 | `hg18 ]
val string_of_genome : genome -> string

val fetch_chrom_size : genome -> string -> unit
val chrom_size_of_file : string -> (string * int) Enum.t










