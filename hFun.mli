type ('a,'b) t

val make : ?init:'a BatEnum.t -> ('a -> 'b) -> ('a, 'b) t
val close : ('a, 'b) t -> ('a, 'b) t

val cache : ('a -> 'b) -> ('a -> 'b)

module Infix : sig
  val ( <-- )  : ('a, 'b) t -> 'a * 'b -> unit
  val ( $ ) : ('a, 'b) t -> 'a -> 'b
end


(* val of_enum : ('a -> 'b) -> 'a Enum.t -> ('a, 'b) t *)
(* val make2 : 'a list -> 'b list -> ('a, 'b) t *)
(* val preimage : ('a -> 'b) -> 'a Enum.t -> ('b, 'a) t *)
(* val compose : ('b -> 'c) -> ('a, 'b) t -> ('a, 'c) t *)
(* val image : ('a, 'b) t -> 'b Enum.t *)
