open Batteries
open Oregon
open Sle.Infix
open Genome

module Neighbouring_peak = struct
  let myaccu () = 
    Accu.create PSet.empty PSet.add

  let shitty_filter s = 
    Enum.fold 
      (fun (seen, accu) (d,r) -> 
	 if PSet.mem r seen then (seen, accu)
	 else (PSet.add r seen, (d,r) :: accu))
      (PSet.empty, [])
      (PSet.enum s)
    |> snd
    |> List.rev

  let run () = 
    let tsss =
      let tss_of_gene = Ucsc.promoter_of_gene `mm9 ~upstream:0 ~downstream:0 in
      let transcripts = Tsv.enum (Ucsc.ensGene `mm9) |> Array.of_enum in
      let t2g = Sle.hrel_of_enum (fun g -> g#transcript, g#gene) (Tsv.enum (Ucsc.ensGtp `mm9)) in
      Array.map (fun t -> tss_of_gene t, List.hd (t#name --> t2g)) transcripts
    and peaks = Chipseq_regions.selected_regions () |> Array.of_list 
    in
    let assoc = Chr_map.leftjoin ~up:(-500000) ~down:500000 
      fst (fun r -> r.Chipseq_regions.loc)
      tsss peaks in
    let assoc = Array.map (fun (t, rs) -> t, Array.(sub rs 0 (Pervasives.min 5 (length rs)))) assoc in
    let assoc = 
      Array.map (fun (t, rs) -> 
		   t, Array.map (fun r -> Location.dist (fst t) r.Chipseq_regions.loc, r) rs) assoc in
    let g2r = 
      let r = myaccu () in
      Array.iter 
	(fun ((_,g), rs) -> Array.iter (Accu.add r g) rs)
	assoc ;
      r
    in
    (Accu.enum g2r) 
    /@ (fun (g,s) -> g, (shitty_filter s |> List.enum |> Enum.take 5 |> Array.of_enum))
    /@ (fun (g,s) -> 
	  Enum.append 
	    (Enum.singleton g) 
	    (Array.enum s /@ (fun (d,r) -> List.enum [ string_of_int d ; Location.to_string r.Chipseq_regions.loc ]) |> Enum.concat) 
	  |> Array.of_enum) 
    /@ Tsv.string_of_line
    |> File.write_lines "rien.tsv"

end
