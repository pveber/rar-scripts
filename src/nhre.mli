type motif = [`direct | `everted | `inverted] * int
type matches = (motif * [`sense | `antisense] * int * float) list

val motifs : motif list
val string_of_motif : motif -> string 


val scan_string : motif list -> float -> string -> matches
