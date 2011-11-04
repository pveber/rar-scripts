open Batteries

type ('a,'b,'c) t
type ('a,'b,'c) accu = ('a,'b,'c) t

val create : ?n:int -> 'c -> ('b -> 'c -> 'c) -> ('a,'b,'c) t
(** [n] is the approximate size of the domain *)

val add : ('a,'b,'c) t -> 'a -> 'b -> unit

val enum : ('a,'b,'c) t -> ('a * 'c) Enum.t

val get  : ('a,'b,'c) t -> 'a -> 'c

module HRel : sig
  type ('a,'b) t = private ('a,'b,'b list) accu
  val of_enum : ('a -> 'b * 'c) -> 'a Enum.t -> ('b, 'c) t
  val enum : ('a,'b) t -> ('a * 'b list) Enum.t

  module Infix : sig
    val ( <+- ) : ('a, 'b) t -> 'a * 'b -> unit
  end
end

(*
    val ( +++ ) : 'a histogram -> 'a -> unit
    val ( +++* ) : 'a histogram -> 'a * int -> unit
*)

module HFun : sig
  type ('a,'b) t
  val make : ('a -> 'b) -> 'a list -> ('a, 'b) t
  val of_enum : ('a -> 'b) -> 'a Enum.t -> ('a, 'b) t
  val make2 : 'a list -> 'b list -> ('a, 'b) t
  val preimage : ('a -> 'b) -> 'a Enum.t -> ('b, 'a) t
  val compose : ('b -> 'c) -> ('a, 'b) t -> ('a, 'c) t
  val image : ('a, 'b) t -> 'b Enum.t

  module Infix : sig
    val ( <-- )  : ('a, 'b) t -> 'a * 'b -> unit
    val ( --> ) : 'a -> ('a, 'b) t -> 'b
  end
end

