open Batteries

module Selection = struct
  type t = (string, ISet.t) PMap.t

  let of_locations e = Location.(
    let accu = Accu.create ISet.empty (fun loc -> ISet.add_range loc.st loc.ed) in
    Enum.iter 
      (fun loc -> Accu.add accu loc.chr loc)
      e ;
    PMap.of_enum (Accu.enum accu)
  )
end
