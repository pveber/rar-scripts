open Batteries
open Oregon
open Target.Infix
open Genome
open RarGenome

let selection_of_bed bed = 
  Tsv.enum bed
  |> Enum.map (fun x -> x#loc)
  |> RarGenome.Selection.of_locations

let jaccard_sim x y = Selection.(
  let xy = float (length (inter x y))
  and x = float (length x)
  and y = float (length y) in
  xy ** 2. /. (x *. y)
)

let symmatrix_init n f = 
  let mat = Array.make_matrix n n 0. in
  for i = 0 to n - 1 do
    for j = 0 to i do 
      let v = f i j in
      mat.(i).(j) <- v ;
      mat.(j).(i) <- v
    done
  done ;
  mat
  
let sample_sim_matrix samples = 
  let selections = Array.map
    (fun s -> selection_of_bed (B.B.Chipseq.macs_peaks s))
    samples in
  symmatrix_init 
    (Array.length samples)
    (fun i j -> jaccard_sim selections.(i) selections.(j))

let sample_hclust_plot samples path = 
  let rp = R.make () 
  and s = sample_sim_matrix samples in
  let dm = Array.map (Array.map (fun x -> 1. -. x)) s in
  let name s = String.quote (B.B.Chipseq.string_of_sample s) in
  R.pdf rp path ;
  R.matrix rp "dm" dm ;
  R.c rp "d <- as.dist(dm)" ;
  R.c rp "plot(hclust(d),labels=c(%s))" (String.concat "," (Array.map name samples |> Array.to_list)) ;
  R.devoff rp ;
  R.close rp
