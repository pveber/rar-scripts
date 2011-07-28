open Batteries

type ('a,'b,'c) t = {
  table : ('a,'c) Hashtbl.t ;
  zero  : 'c ;
  add   : 'b -> 'c -> 'c
}

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
