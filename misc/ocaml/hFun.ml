type ('a, 'b) t = ('a,'b) Hashtbl.t

module Infix = struct
  let ( <-- ) hfun (x, v) = Hashtbl.replace hfun x v
    
  let ( $ ) hfun x = Hashtbl.find hfun x
end

open Infix

let make ?(init = BatEnum.empty ()) f = 
  let r  = Hashtbl.create 251 in
  BatEnum.iter 
    (fun x -> r <-- (x, f x))
    init ;
  r

let of_enum = BatHashtbl.of_enum

let domain hfun = 
  BatHashtbl.keys hfun

let image hfun = 
  BatHashtbl.values hfun

let cache f = 
  let hfun = make f in
  fun x -> 
    try Hashtbl.find hfun x
    with Not_found -> (
      let y = f x in
      Hashtbl.replace hfun x y ;
      y
    )

(* let of_enum f support =  *)
(*   let r = Hashtbl.create 1024 in *)
(*   Enum.iter (fun x -> r <-- (x, f x)) support ; *)
(*   r *)

(* let make2 domain image =  *)
(*   let n = List.length domain  *)
(*   and m = List.length image in *)
(*   let _ = if n <> m then raise (Invalid_argument "hfun_extension") in *)
(*   let r = Hashtbl.create n in *)
(*   List.iter2 (fun x y -> r <-- (x, y)) domain image ; *)
(*   r *)

(* let preimage f support =  *)
(*   let r = Hashtbl.create 251 in *)
(*   Enum.iter (fun x -> r <-- (f x, x)) support ; *)
(*   r *)

(* let compose f hfun = Hashtbl.map (fun _ v -> f v) hfun *)

(* let image hfun = Hashtbl.values hfun *)
