open Batteries

val shell : string -> unit
(** Shell invocation. Raises failwith in case of non zero return code *)

val sh : ('a,unit,string,unit) format4 -> 'a
(** Same as [shell], printf-style *)

val pipefail : string -> string -> unit
(** [pipefail cmd1 cmd2] does a special pipe of commands [cmd1] and [cmd2]
    whose exit code is non null as soon as that of any of the two commands 
    is. 
    If the exist code is non zero, calls [failwith]
*)

val save : string -> 'a -> unit
val load : string -> 'a

val with_tmp_filename : (string -> unit) -> unit

val string_split : char -> string -> string array
val string_split_noeps : char -> string -> string array
val string_quote : string -> string
val string_unquote : string -> string
