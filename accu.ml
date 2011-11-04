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



module HFun = struct
  type ('a,'b) t = ('a,'b) Hashtbl.t

  module Infix = struct
    let ( <-- ) f (x, v) = Hashtbl.replace f x v
	   
    let ( --> ) x f = Hashtbl.find f x
  end

  open Infix

  let make f support = 
    let r = Hashtbl.create (List.length support) in
    List.iter (fun x -> r <-- (x, f x)) support ;
    r

  let of_enum f support = 
    let r = Hashtbl.create 1024 in
    Enum.iter (fun x -> r <-- (x, f x)) support ;
    r

  let make2 domain image = 
    let n = List.length domain 
    and m = List.length image in
    let _ = if n <> m then raise (Invalid_argument "hfun_extension") in
    let r = Hashtbl.create n in
    List.iter2 (fun x y -> r <-- (x, y)) domain image ;
    r

  let preimage f support = 
    let r = Hashtbl.create 251 in
    Enum.iter (fun x -> r <-- (f x, x)) support ;
    r

  let compose f hfun = Hashtbl.map (fun _ v -> f v) hfun

  let image hfun = Hashtbl.values hfun
end

(*
    let ( +++ ) c k =
      let n = try Hashtbl.find c k with Not_found -> 0 in
      Hashtbl.replace c k (n + 1)

    let ( +++* ) c (k,i) =
      let n = try Hashtbl.find c k with Not_found -> 0 in
      Hashtbl.replace c k (n + i)
*)
