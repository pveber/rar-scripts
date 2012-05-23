(** A datatype to cache intermediary results on files *)

type 'a t
type 'a serializer = (string -> 'a -> unit) * (string -> 'a)

val value : ?serializer:'a serializer -> string -> (unit -> 'a) -> 'a t
(** [value fn f] caches the result of [f ()] in a file at [fn]. The
    default serializer is from the [Marshal] module *)


val file : string -> (string -> unit) -> string t

val get : 'a t -> 'a





















