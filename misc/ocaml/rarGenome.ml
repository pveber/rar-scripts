open Batteries
open Genome



(*
module Location = struct
  open Batteries
  open Printf

  type t = {
    chr : string ;
    st : int ;
    ed : int
  }

  let to_string l = sprintf "%s:%d-%d" l.chr l.st l.ed

  let length l = l.ed - l.st + 1

  let make chr st ed = 
    if st > ed then (
      let msg = sprintf "Genome.Location.make: incorrect coordinates %d %d\n" st ed in
      raise (Invalid_argument msg) ;
    ) ;
    { chr = chr ; st = st ; ed = ed }

  let of_string s = 
    Scanf.sscanf 
      s "%s@:%d@-%d" 
      make

  let upstream ~up ~down strand loc = match strand with
    | `Sense ->  
      make loc.chr (loc.st - up + 1) (loc.st + down)
    | `Antisense -> 
      make loc.chr (loc.ed - down + 1) (loc.ed + up)

  let relmove l st ed = { l with st = l.st + st ; ed = l.ed + ed }

  let included_in s s' =
    if s.chr <> s'.chr then false 
    else (s.st >= s'.st && s.ed <= s'.ed)

  let intersection s s' = 
    if s.chr <> s'.chr then false else (
      let p s s' = (s.ed >= s'.st) && (s.ed <= s'.ed)
      in p s s' || p s' s
    )

  let inter s s' = 
    if not (intersection s s') then raise (Invalid_argument "Genome.Location.inter")
    else (make s.chr (max s.st s'.st) (min s.ed s'.ed))
      
  let dist s s' = 
    if s.chr <> s'.chr then raise (Invalid_argument "Genome.Location.dist")
    else (
      if intersection s s' then 0
      else min (abs (s'.st - s.ed)) (abs (s.st - s'.ed))
    )

  let position ~from loc = 
    if loc.chr <> from.chr then raise (Invalid_argument "Genome.Location.position")
    else (
      if intersection from loc then 0
      else (
	let a, b = from.st - loc.ed, loc.st - from.ed in
	if abs a < abs b then a else b
      )
    )
      
  let compare (x : t) (y : t) = compare x y    
end
*)

let location_upstream ~up ~down strand loc = Location.(
  match strand with
    | `Sense ->  
	make loc.chr (max 0 (loc.st - up)) (max 0 (loc.st + down))
    | `Antisense -> 
	make loc.chr (max 0 (loc.ed - down)) (max 0 (loc.ed + up))
)

module Selection = struct
  module M = Map.StringMap
  type t = ISet.t M.t

  let of_locations e = Location.(
    let accu = Accu.create ISet.empty (fun loc -> ISet.add_range loc.st loc.ed) in
    Enum.iter (fun loc -> Accu.add accu loc.chr loc) e ;
    M.of_enum (Accu.enum accu)
  )

  let inter u v =
    M.fold
      (fun k set_u accu ->
	 try 
	   let set_v = M.find k v in
	   M.add k (ISet.inter set_u set_v) accu
	 with Not_found -> accu)
      u M.empty

  let diff u v =
    M.fold
      (fun k set_u accu ->
	 let set_u' = 
	   try 
	     let set_v = M.find k v in
	     ISet.diff set_u set_v
	   with Not_found -> set_u
	 in M.add k set_u' accu)
      u M.empty

  let length x = 
    M.fold (fun _ set accu -> ISet.cardinal set + accu) x 0

  let qinclusion sel loc = Location.(ISet.(
    try 
      let k = 
	inter
	  (add_range loc.st loc.ed empty)
	  (M.find loc.chr sel)
        |> cardinal in
      float k /. float (loc.ed - loc.st + 1)
    with Not_found -> 0.
  ))
end



