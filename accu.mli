open Batteries

type ('a,'b,'c) t

val create : ?n:int -> 'c -> ('b -> 'c -> 'c) -> ('a,'b,'c) t
(** [n] is the approximate size of the domain *)

val add : ('a,'b,'c) t -> 'a -> 'b -> unit

val enum : ('a,'b,'c) t -> ('a * 'c) Enum.t

val get  : ('a,'b,'c) t -> 'a -> 'c
