open Batteries

module Selection : sig
  type t = (string, ISet.t) PMap.t
  val of_locations : Location.t Enum.t -> t
end
