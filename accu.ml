open Batteries

type ('a,'b,'c) t = {
  table : ('a,'c) Hashtbl.t ;
  zero  : 'c ;
  add   : 'b -> 'c -> 'c
}

type ('a,'b,'c) accu = ('a,'b,'c) t

let create ?(n = 251) zero f = {
  table = Hashtbl.create 251 ;
  zero = zero ;
  add = f
}

let add t x y = 
  let accu = 
    try Hashtbl.find t.table x
    with Not_found -> t.zero in
  Hashtbl.replace t.table x (t.add y accu)

let enum t = Hashtbl.enum t.table

let get t = Hashtbl.find t.table


module HRel = struct
  type ('a,'b) t = ('a,'b,'b list) accu

  module Infix = struct
    let ( <+- ) r (k, v) = add r k v
  end

  open Infix

  let of_enum f support = 
    let r = create [] (fun h t -> h :: t) in
    BatEnum.iter (fun z -> let x,y = f z in add r x y) support ;
    r

  let enum = enum
end



let cache f = 
  let t = Hashtbl.create 251 in 
  fun x -> 
    try Hashtbl.find t x 
    with Not_found -> (
      let y = f x in 
      Hashtbl.add t x y ;
      y
    )

(*
    let ( +++ ) c k =
      let n = try Hashtbl.find c k with Not_found -> 0 in
      Hashtbl.replace c k (n + 1)

    let ( +++* ) c (k,i) =
      let n = try Hashtbl.find c k with Not_found -> 0 in
      Hashtbl.replace c k (n + i)
*)
