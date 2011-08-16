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
  
let sim_matrix tracks = 
  symmatrix_init 
    (Array.length tracks)
    (fun i j -> jaccard_sim (snd tracks.(i)) (snd tracks.(j)))

let dist_of_sim m = 
  Array.map (Array.map (fun x -> 1. -. x)) m

let dist_matrix = sim_matrix |- dist_of_sim

let tracks_of_samples samples = 
  Array.map
    (fun s -> 
       B.B.Chipseq.(string_of_sample s,
		    (macs_peaks |- Macs.filter ~pvalue:1e-9 |- selection_of_bed) s))
    samples

type col es_tsv = {
  loc : Location ;
  x : int
}

let es_peaks path = 
  Misc.wget 
    ~gunzip:true
    ("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM288nnn/" ^ path)
    Tsv.({has_header = false ; parse = Es_tsv.of_array })

let oct4_peaks = es_peaks "GSM288346/GSM288346_ES_Oct4.txt.gz"

let sample_hclust_plot tracks path = 
  let rp = R.make () 
  and dm = dist_matrix tracks
  and labels = 
    Array.map
      (fst |- String.quote)
      tracks
    |> Array.to_list
    |> String.concat "," in
  R.pdf rp path ;
  R.matrix rp "dm" dm ;
  R.c rp "d <- as.dist(dm)" ;
  R.c rp "plot(hclust(d),labels=c(%s))" labels ;
  R.devoff rp ;
  R.close rp
