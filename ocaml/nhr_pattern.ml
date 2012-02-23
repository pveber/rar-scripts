open Batteries
open Printf
open Biocaml

module A = Archimedes

let fa_input = Sys.argv.(1)
let output = Sys.argv.(2)
let _ = ignore (Sys.command ("mkdir -p " ^ output))
let theta = Sys.argv.(3)

let selected_motifs = List.concat [
  List.init 10 (fun i -> `direct, i) ;
  List.init 10 (fun i -> `everted, i) ;
  List.init 10 (fun i -> `inverted, i) ;
]


(* Matrix-related functions *)
let margin f init mat = 
  Array.fold_left 
    (Array.fold_left f)
    init mat

let cmargin f init mat =
  let n = Array.length mat.(0) in
  let r = Array.make n init in
  for j = 0 to n - 1 do
    for i = 0 to Array.length mat - 1 do
      r.(j) <- f r.(j) mat.(i).(j)
    done
  done ;
  r

let init_matrix n m f = Array.init n (fun i -> Array.init m (fun j -> f i j))

let mat_map2 f mat1 mat2 = 
  let m = Array.length mat1 
  and n = Array.length mat2 in
  init_matrix m n (fun i j -> f mat1.(i).(j) mat2.(i).(j))

let sum  = Array.fold_left (+)  0
let sumf = Array.fold_left (+.) 0.

(* Statistics *)
let function_of_enum neutral e = 
  let c = Hashtbl.of_enum e in
  fun x -> Hashtbl.find_default c x neutral

let motif_occurrence_vector occ =
  List.enum occ
  |> Biocaml.Accu.counts Tuple5.first
  |> function_of_enum 0
  

let motif_counts_byseq occbyseq = 
  let motifs = Array.of_list selected_motifs in
  Array.map
    (fun occ -> Array.map (motif_occurrence_vector occ) motifs)
    occbyseq

let motif_counts occbyseq = 
  cmargin (+) 0 (motif_counts_byseq occbyseq)

let related_motifs ((_,_,st,ed,_),(_,_,st',ed',_)) = 
  (ed = st' + 6) || (ed' = st + 6) || 
  (st = st' && ed <> ed') || (ed = ed' && st <> st')

let inclusion_link ((_,_,st,ed,_),(_,_,st',ed',_)) = 
  (st = st' && ed <> ed') || (ed = ed' && st <> st')

let proximity_link ((_,_,st,ed,_),(_,_,st',ed',_)) = 
  (ed = st' + 6) || (ed' = st + 6)

let motif_motif_mat_cm_for_a_seq f motifs occ =
  List.cartesian_product occ occ
  |> List.filter f
  |> List.enum
  |> Biocaml.Accu.counts (fun ((x,_,_,_,_),(y,_,_,_,_)) -> x,y)
  |> function_of_enum 0  

let motif_motif_cm_byseq f occbyseq = 
  let motifs = Array.of_list selected_motifs in
  Array.map 
    (fun occ -> 
       let mat = motif_motif_mat_cm_for_a_seq f motifs occ in
       let n = Array.length motifs 
       and f i j = mat (motifs.(i),motifs.(j)) in
       init_matrix n n f)
    occbyseq

let motif_motif_mat_cm f occbyseq = 
  Array.fold_left
    (mat_map2 (+))
    (let n = List.length selected_motifs in Array.make_matrix n n 0)
    (motif_motif_cm_byseq f occbyseq)

let double f x = f x x

let mi nij ni nj n = 
  let pij = float nij /. float n
  and pi  = float ni  /. float n
  and pj  = float nj  /. float n in
  pij *. log (pij /. (pi *. pj))

open Graph

module V = struct
  type t = int
  let compare = compare
  let equal = ( = )
  let hash = Hashtbl.hash
end

module E = struct
  type kind = Inclusion | Proximity
  type t = kind * float
  let compare = compare
  let default = Proximity, 0.
end

module G = Persistent.Graph.ConcreteLabeled(V)(E)
module Op = Oper.P(G)

let association_graph theta kind f occbyseq =
  let motifs = Array.of_list selected_motifs in
  let pair_counts = motif_motif_mat_cm f occbyseq in
  let margins = Array.map sum pair_counts in
  let total = sum margins in
  List.(Array.range motifs |> of_enum |> double cartesian_product |> enum)
  // (fun (i,j) -> i >= j)
  /@ (fun (i,j) -> mi pair_counts.(i).(j) margins.(i) margins.(j) total, i, j)
  // (fun (x,_,_) -> x > theta)
  |> Enum.fold (fun accu (v,x,y) -> G.add_edge_e accu G.E.(create x (kind,v) y)) G.empty
      
let global_association_graph occbyseq =
  let g_incl = association_graph 0.01 E.Inclusion inclusion_link occbyseq
  and g_prox = association_graph 0.004 E.Proximity proximity_link occbyseq in
  Op.union g_incl g_prox


    
let association_graph_neato_output motif_freq g fn = 
  let motifs = Array.of_list selected_motifs in
  let module G = struct
    let graph_attributes g = []
    let default_vertex_attributes g = []
    let vertex_name i = Nhre.string_of_motif motifs.(i)
    let vertex_attributes i = [ `Shape `Circle ; `Height (motif_freq.(i) ** 3.) ; ]
    let get_subgraph _ = None
    let default_edge_attributes _ = [ `Len 2. ]
    let edge_attributes e = 
      let kind, weight = G.E.label e in [
	`Color (match kind with E.Inclusion -> 0xff0000 | E.Proximity -> 0x00ff00) ;
	`Style (if weight > 0.05 then `Bold else if weight > 0.01 then `Solid else `Dashed)
      ]
    include G
  end in
  let module N = Graphviz.Neato(G) in
  let oc = Pervasives.open_out fn in
  N.output_graph oc g ;
  Pervasives.close_out oc


let association_graph_neato_layout motif_freq g =
  let module LayoutNode = struct
    type t = int * (float * float)
    let compare = compare
    let equal = ( = )
    let hash = Hashtbl.hash
  end in
  let module LayoutGraph = Persistent.Graph.Concrete(LayoutNode) in
  let module LabelBuilder = struct
    open Dot_ast
    let motif_ids = Array.(map Nhre.string_of_motif (of_list selected_motifs))

    let soi = function
      | Ident s -> s
      | Html s -> s
      | Number s -> s
      | String s -> s

    let aoa attrs = 
      List.map
	(List.map (fun (k,v) -> soi k, Option.map_default soi "" v))
	attrs
      |> List.concat

    let pos s = Scanf.sscanf s "%f,%f" (fun x y -> x,y)

    let node (id,_) attrs =
      (Array.findi (( = ) (soi id)) motif_ids,
       pos (List.assoc "pos" (aoa attrs)))

    let edge attrs = ()
  end in
  let module P = Dot.Parse(Builder.P(LayoutGraph))(LabelBuilder) in
  association_graph_neato_output motif_freq g "delme.neato" ;
  ignore (Sys.command "neato delme.neato > delme2.neato") ;
  let g,bbox,_ = P.parse_bounding_box_and_clusters "delme2.neato" in
  LayoutGraph.fold_vertex
    (fun h t -> h :: t)
    g [],
  Scanf.sscanf bbox "%f,%f,%f,%f" (fun _ _ x y -> x, y)
  
let pi = 4. *. atan 1.

module GraphDraw = struct
  open Archimedes
  open Viewport

  let circle ?(filled = false) vp x y r = 
    move_to vp (x +. r) y ;
    arc vp r 0. (2. *. pi) ;
    (if filled then fill else stroke) vp `Device

  let node ?(r = 0.02) ?(bg = Color.white) ?(fg = Color.black) vp x y label = 
    set_line_width vp 1. ;
    set_color vp bg ;
    circle ~filled:true vp x y r ;
    set_color vp fg ;
    circle vp x y r ;
    text ~coord:`Device vp x y label

  let edge ?(width = 1.) ?(col = Color.black) vp x y x' y' =
    set_color vp col ;
    set_line_width vp width ;
    move_to vp x y ;
    line_to vp x' y' ;
    stroke vp `Device
end

let association_graph_archimedes_output motif_freq (w,h) node_layout g = 
  let h = max w h
  and w = max w h in
  let vp = A.init ~w ~h (* ["graphics" ; "hold"]*) ["cairo"; "PDF"; "rien.pdf"] in
  GraphDraw.(
    G.iter_edges_e
      (fun (u,l,v) -> 
	let x_u, y_u = List.assoc u node_layout 
	and x_v, y_v = List.assoc v node_layout 
	and col, width = match l with 
	    E.Inclusion, weight -> A.Color.red, 1. +. 100. *. weight
	  | E.Proximity, weight -> A.Color.blue, 1. +. 100. *. weight
	in
	edge vp ~col ~width (x_u /. w) (y_u /. h) (x_v /. w) (y_v /. h))
      g ;
    G.iter_vertex 
      (fun v -> 
	let x,y = List.assoc v node_layout 
	and r = 0.02 +. motif_freq.(v) /. 15. in
	node vp ~r (x /. w) (y /. h) (Nhre.string_of_motif (List.nth selected_motifs v)))
      g
  ) ;
  A.close vp

  

let association_report label f occbyseq = 
  let motifs = Array.of_list selected_motifs in
  let pair_counts = motif_motif_mat_cm f occbyseq in
  let margins = Array.map sum pair_counts in
  let total = sum margins in
  printf "%s criterion\n" label ;
  List.(Array.range motifs |> of_enum |> double cartesian_product |> enum)
  // (fun (i,j) -> i >= j)
  /@ (fun (i,j) -> mi pair_counts.(i).(j) margins.(i) margins.(j) total, motifs.(i), motifs.(j))
  |> List.of_enum
  |> List.sort (fun x y -> compare y x)
  |> List.take 30
  |> List.iter (fun (score,mot_i,mot_j) -> 
		  printf "  %g\t%s\t%s\n" score (Nhre.string_of_motif mot_i) (Nhre.string_of_motif mot_j))


let fa =
  Fasta.enum_of_file fa_input
  |> snd
  |> List.of_enum
(*  |> List.take 100 *)

let sequences = List.map snd fa

let motif_search theta seq = 
  Nhre.matches_of_sequence selected_motifs Nhre.estimated_counts theta seq
  |> List.map (fun (m,sense,pos,score) -> (m,sense,pos,pos + Nhre.motif_length m,score))

let occbyseq theta = Array.(map (motif_search theta) (of_list sequences))
  

let () = 
  let occbyseq = occbyseq (float_of_string theta) in
  let n = Array.length occbyseq in
  let motif_counts_byseq = motif_counts_byseq occbyseq in
  let motif_counts = cmargin (+) 0 motif_counts_byseq in
  let motif_freq = Array.map (fun c -> float c /. float n) motif_counts in
  printf "Number of sequences: %d\n" n ;
  printf "Total number of matches: %d\n" (margin (+) 0 motif_counts_byseq) ;
  printf "Motif frequency:\n" ;
  List.iteri
    (fun i motif ->
       printf "  %s : %g\n" (Nhre.string_of_motif motif) motif_freq.(i))
    selected_motifs ;
  printf "Association strength:\n" ;
  association_report "Hexamer-sharing" related_motifs occbyseq ;
  association_report "Inclusion link" inclusion_link occbyseq ;
  association_report "Proximity link" proximity_link occbyseq ;
  let g = global_association_graph occbyseq in
  let node_layout,bbox = association_graph_neato_layout motif_freq g in
  association_graph_archimedes_output motif_freq bbox node_layout g
