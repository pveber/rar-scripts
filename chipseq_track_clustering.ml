open Batteries
open Printf
open Oregon
open Target.Infix
open Genome
open RarGenome

let sh = Sle.sh

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

let recenter n l = Location.(
  let i = (l.st + l.ed) / 2 in
  move l (i - n) (i + n)
)

let macs_recenter n p = Location.(
  let i = p#loc.st + p#summit in
  move p#loc (i - n) (i + n)
)

let tracks_of_samples radius samples = 
  Array.map
    (fun s -> 
       B.B.Chipseq.(string_of_sample s,
		    (macs_peaks |- Macs.filter ~pvalue:1e-9 |- Macs.locations_around_summit ~radius |- selection_of_bed) s))
    (Array.of_list samples)

type col es_tsv = {
  loc : Location ;
  x : int
}

let es_peaks path = 
  Misc.wget 
    ~gunzip:true
    ("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM288nnn/" ^ path)
    Tsv.({has_header = false ; parse = Es_tsv.of_array })

let nanog_peaks = es_peaks "GSM288345/GSM288345_ES_Nanog.txt.gz"
let oct4_peaks = es_peaks "GSM288346/GSM288346_ES_Oct4.txt.gz"
let sox2_peaks = es_peaks "GSM288347/GSM288347_ES_Sox2.txt.gz"
let smad1_peaks = es_peaks "GSM288348/GSM288348_ES_Smad1.txt.gz"
let e2f1_peaks = es_peaks "GSM288349/GSM288349_ES_E2f1.txt.gz"
let tcfcp2l1_peaks = es_peaks "GSM288350/GSM288350_ES_Tcfcp2l1.txt.gz"
let zfx_peaks = es_peaks "GSM288352/GSM288352_ES_Zfx.txt.gz"
let stat3_peaks = es_peaks "GSM288353/GSM288353_ES_Stat3.txt.gz"
let klf4_peaks = es_peaks "GSM288354/GSM288354_ES_Klf4.txt.gz"
let c_myc_peaks = es_peaks "GSM288356/GSM288356_ES_c-Myc.txt.gz"
let n_myc_peaks = es_peaks "GSM288357/GSM288357_ES_n-Myc.txt.gz"
let p300_peaks = es_peaks "GSM288359/GSM288359_ES_p300.txt.gz"
let suz12_peaks = es_peaks "GSM288360/GSM288360_ES_Suz12.txt.gz"


let es_selection n x = 
  Tsv.enum x 
  |> Enum.map (fun x -> recenter n x.loc)
  |> RarGenome.Selection.of_locations

let es_tracks n = [|
  "Nanog", es_selection n nanog_peaks ;
  "Oct4", es_selection n oct4_peaks ;
  "Sox2", es_selection n sox2_peaks ;
  "Smad1", es_selection n smad1_peaks ;
  "E2f1", es_selection n e2f1_peaks ;
  "Tcfcp2l1", es_selection n tcfcp2l1_peaks ;
  "Zfx", es_selection n zfx_peaks ;
  "Stat3", es_selection n stat3_peaks ;
  "Klf4", es_selection n klf4_peaks ;
  "c-Myc", es_selection n c_myc_peaks ;
  "n-Myc", es_selection n n_myc_peaks ;
  "p300", es_selection n p300_peaks ;
  "Suz12", es_selection n suz12_peaks ;
|]

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
