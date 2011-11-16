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
  xy *.2. /. (x +. y)
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
  move l (max 0 (i - n)) (i + n)
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


module Wei_dataset = struct

  type col tsv_format = {
    loc : Location ;
    x : int
  }

  let peaks path = 
    Misc.wget 
      ~gunzip:true
      ("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM288nnn/" ^ path)
      Tsv.({has_header = false ; parse = Tsv_format.of_array })

  let nanog_peaks = peaks "GSM288345/GSM288345_ES_Nanog.txt.gz"
  let oct4_peaks = peaks "GSM288346/GSM288346_ES_Oct4.txt.gz"
  let sox2_peaks = peaks "GSM288347/GSM288347_ES_Sox2.txt.gz"
  let smad1_peaks = peaks "GSM288348/GSM288348_ES_Smad1.txt.gz"
  let e2f1_peaks = peaks "GSM288349/GSM288349_ES_E2f1.txt.gz"
  let tcfcp2l1_peaks = peaks "GSM288350/GSM288350_ES_Tcfcp2l1.txt.gz"
  let zfx_peaks = peaks "GSM288352/GSM288352_ES_Zfx.txt.gz"
  let stat3_peaks = peaks "GSM288353/GSM288353_ES_Stat3.txt.gz"
  let klf4_peaks = peaks "GSM288354/GSM288354_ES_Klf4.txt.gz"
  let c_myc_peaks = peaks "GSM288356/GSM288356_ES_c-Myc.txt.gz"
  let n_myc_peaks = peaks "GSM288357/GSM288357_ES_n-Myc.txt.gz"
  let p300_peaks = peaks "GSM288359/GSM288359_ES_p300.txt.gz"
  let suz12_peaks = peaks "GSM288360/GSM288360_ES_Suz12.txt.gz"


  let selection n x = 
    Tsv.enum x 
    |> Enum.map (fun x -> recenter n x.loc)
    |> RarGenome.Selection.of_locations

  let tracks n = [|
    "Wei Nanog", selection n nanog_peaks ;
    "Wei Oct4", selection n oct4_peaks ;
    "Wei Sox2", selection n sox2_peaks ;
    "Wei Smad1", selection n smad1_peaks ;
    "Wei E2f1", selection n e2f1_peaks ;
    "Wei Tcfcp2l1", selection n tcfcp2l1_peaks ;
    "Wei Zfx", selection n zfx_peaks ;
    "Wei Stat3", selection n stat3_peaks ;
    "Wei Klf4", selection n klf4_peaks ;
    "Wei c-Myc", selection n c_myc_peaks ;
    "Wei n-Myc", selection n n_myc_peaks ;
    "Wei p300", selection n p300_peaks ;
    "Wei Suz12", selection n suz12_peaks ;
  |]
end
  
let rar_es url = Target.F.make 
  (object
     method id = sprintf "Chipseq_track_clustering.rar_es[%s]" url
     method deps = []
     method build path = 
       sh "mkdir -p %s.tmp" path ;
       sh "wget -O %s.tmp/GSM.gz %s" path url ;
       sh "gunzip %s.tmp/GSM.gz" path ;
       sh "mv %s.tmp/GSM %s" path path ;
       sh "rm -rf %s.tmp" path
     method ty = Tsv.({ has_header = true ; parse = (fun x -> Location.of_string ("chr" ^ x.(0))) })
   end)

let rar_es_d2 = rar_es "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM482nnn/GSM482749/GSM482749%5FRARd2%5Fvs%5FWCE%2Epeaks%2Etxt%2Egz"

let rar_es_d2_8 = rar_es "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM482nnn/GSM482750/GSM482750%5FRARd2p8h%5Fvs%5FWCE%2Epeaks%2Etxt%2Egz"

let rar_es_selection n x = 
  Tsv.enum x 
  |> Enum.map (fun x -> recenter n x)
  |> RarGenome.Selection.of_locations

let rar_es_tracks n = [|
   "RAR-ES D2", rar_es_selection n rar_es_d2 ;
   "RAR-ES D2 + 8", rar_es_selection n rar_es_d2_8 ;
|]



let smad2_es url = Target.F.make 
  (object
     method id = sprintf "Chipseq_track_clustering.smad2_es[%s]" url
     method deps = []
     method build path = 
       sh "mkdir -p %s.tmp" path ;
       sh "wget -O %s.tmp/GSM.gz %s" path url ;
       sh "gunzip %s.tmp/GSM.gz" path ;
       sh "mv %s.tmp/GSM %s" path path ;
       sh "rm -rf %s.tmp" path
     method ty = Tsv.({ has_header = true ; parse = (fun x -> Location.make x.(0) (int_of_string x.(1)) (int_of_string x.(1))) })
   end)

let smad2_es_18h_activin = smad2_es "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM578nnn/GSM578474/GSM578474%5FACT%5FIPvsWCE%5Fpeaks%2Etxt%2Egz"
let smad2_es_18h_sb = smad2_es "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM578nnn/GSM578476/GSM578476%5FSB%5FIPvsWCE%5Fpeaks%2Etxt%2Egz"
let smad2_es_18h_dmso = smad2_es "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM578nnn/GSM578475/GSM578475%5FDMSO%5FIPvsWCE%5Fpeaks%2Etxt%2Egz"

let smad2_es_selection n x =
  Tsv.enum x
  |> Enum.map (fun x -> recenter n x)
  |> RarGenome.Selection.of_locations

let smad2_es_tracks n = [|
   "pSMAD2-ES 18h + activin", smad2_es_selection n smad2_es_18h_activin ;
   "pSMAD2-ES 18h + SB", smad2_es_selection n smad2_es_18h_sb ;
   "pSMAD2-ES 18h + DMSO", smad2_es_selection n smad2_es_18h_dmso ;
|]

let gsm_fetch ty url = Target.F.make
  (object
     method id = sprintf "Chipseq_track_clustering.gsm_fetch[%s]" url
     method deps = []
     method build path =
       sh "mkdir -p %s.tmp" path ;
       sh "wget -O %s.tmp/GSM.gz %s" path url ;
       sh "gunzip %s.tmp/GSM.gz" path ;
       sh "mv %s.tmp/GSM %s" path path ;
       sh "rm -rf %s.tmp" path
     method ty = ty
   end)


module Schnetz_dataset = struct
  (* http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22341 *)

  let ty = Tsv.({ has_header = false ; 
		  parse = new Bed.base_parser })

  let p300_high = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558675/GSM558675%5Fp300%5Fpeak%5FHigh%5Fthresh%2Etxt%2Egz"
  let p300_low = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558675/GSM558675%5Fp300%5Fpeak%5FLow%5Fthresh%2Etxt%2Egz"
  let p300_middle = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558675/GSM558675%5Fp300%5Fpeak%5FMiddle%5Fthresh%2Etxt%2Egz"

  let chd7_high = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558674/GSM558674%5FCHD7%5Fpeak%5FHigh%5Fthresh%2Etxt%2Egz"
  let chd7_low = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558674/GSM558674%5FCHD7%5Fpeak%5FLow%5Fthresh%2Etxt%2Egz"
  let chd7_middle = gsm_fetch ty "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM558nnn/GSM558674/GSM558674%5FCHD7%5Fpeak%5FMiddle%5Fthresh%2Etxt%2Egz"

  let selection n x =
    Tsv.enum x
    |> Enum.map (fun x -> recenter n x#loc)
    |> RarGenome.Selection.of_locations

  let tracks n = [|
    "Schnetz p300 high", selection n p300_high ;
    "Schnetz p300 low", selection n p300_low ;
    "Schnetz p300 middle", selection n p300_middle ;
    "Schnetz Chd7 high", selection n chd7_high ;
    "Schnetz Chd7 low", selection n chd7_low ;
    "Schnetz Chd7 middle", selection n chd7_middle ;
  |]

end

let encode_tsv_parser l = 
  Location.make l.(1) (int_of_string l.(2)) (int_of_string l.(3))

let p300_encode_tsv =
  Target.F.input "manual/chipseq/peaks/p300_encode.tsv" Tsv.({ has_header = false ; parse = encode_tsv_parser })

let p300_encode n = 
  ("p300 ENCODE",
   rar_es_selection n p300_encode_tsv)

let all_tracks n = Array.concat [
  Wei_dataset.tracks n ;
  tracks_of_samples n Selected_chipseq_regions.conditions ;
  rar_es_tracks n ;
  smad2_es_tracks n ;
  Schnetz_dataset.tracks n ;
  [| p300_encode n |]
]

let save_dist_matrix tracks path = 
  let rp = R.make () 
  and dm = dist_matrix tracks
  and labels = 
    Array.map
      (fst |- String.quote)
      tracks
    |> Array.to_list
    |> String.concat "," in
  R.matrix rp "dm" dm ;
  R.c rp "d <- as.dist(dm)" ;
  R.c rp "labels <- c(%s)" labels ;
  R.c rp "r <- list(dist = d, labels = labels)" ;
  R.c rp "save(r,file = '%s')" path ;
  R.close rp

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


(* works for Location.t Tsv.file *)
let bed_of_dataset x = Target.F.make
  (object
     method id = "Chipseq_track_clustering.bed_of_dataset[r1]"
     method deps = [] ++ x
     method build path = Tsv.enum x
       |> Enum.map (Bed.unparser |- Tsv.string_of_line)
       |> File.write_lines path
     method ty = Tsv.({ has_header = false ; parse = new Bed.base_parser })
   end)
  