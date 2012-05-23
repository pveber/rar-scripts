open Batteries

val shell : string -> unit
(** Shell invocation. Raises failwith in case of non zero return code *)

val sh : ('a,unit,string,unit) format4 -> 'a
(** Same as [shell], printf-style *)

val save : string -> 'a -> unit
val load : string -> 'a

val with_tmp_filename : (string -> unit) -> unit




















