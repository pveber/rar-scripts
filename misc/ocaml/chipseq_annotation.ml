open Batteries
open Oregon
open Genome
open Oregon.Target
open Oregon.Target.Infix
open HFun.Infix

type annot = {
  loc : Location.t ;
  counts : int array ;
  coverage : float array ; 
  pvalue : float array ;
}
type t = annot array * B.Chipseq.sample array * (B.Chipseq.sample * B.Chipseq.sample) array

let pmt_test x1 x2 r1 r2 = 
  let l1 = float x1 /. float r1 
  and l2 = float x2 /. float r2 in 
  if l1 < l2 then Pmt.log_poisson_margin x1 x2 r1 r2
  else Pmt.log_poisson_margin x2 x1 r2 r1

let make design locations = 
  let all_samples = 
    List.fold_left (fun accu (x,y) -> x :: y :: accu) [] design
    |> List.sort_unique compare 
  and bowtie = HFun.cache B.Chipseq.bowtie_wodup in
  let chrmap = HFun.cache (bowtie |- Sam.to_chr_map)
  and nbmappings = HFun.cache (bowtie |- Sam.nbmappings) in
  V.make (object
    method id = "Chipseq_annotation[r4]" 
    method deps = [] ++ locations ++* (List.map chrmap all_samples) ++* (List.map nbmappings all_samples)
    method build = 
      let locations = locations#value in 
      let _ = Pmt.prime_cache() in
      let counts_and_coverage_by_sample = 
        List.map 
          (fun x -> 
            Gc.major () ; Gc.major () ;
            let read_map = (chrmap x)#value in
            let library_size = (nbmappings x)#value in
            let counts = Array.map (Chr_map.nbregions_in read_map) locations in
            let coverage =
              Array.mapi
                (fun i loc ->
                  float counts.(i) /. float library_size *. 1e6 /. (loc |> Location.length |> float) *. 1000.)
                locations
            in x, (counts, coverage))
          all_samples 
      in
      let counts = HFun.cache (fun x -> fst (List.assoc x counts_and_coverage_by_sample))
      and coverage = HFun.cache (fun x -> snd (List.assoc x counts_and_coverage_by_sample)) 
      and library_size = HFun.cache (fun x -> (nbmappings x)#value)
      and all_samples = Array.of_list all_samples
      and design = Array.of_list design
      in
      Array.mapi
        (fun i loc ->
          let counts = Array.map   (fun x -> (counts x).(i)) all_samples
          and coverage = Array.map (fun x -> (coverage x).(i)) all_samples
          and pvalue = Array.map (fun (x,y) -> -. pmt_test (counts x).(i) (counts y).(i) (library_size x) (library_size y)) design
          in
          { loc ; counts ; coverage ; pvalue })
        locations,
      all_samples,
      design
  end)

let tsv output_path (annot,domain,domain2) = 
  let labels l = 
    Array.map (fun x -> l ^ " " ^ (B.Chipseq.string_of_sample x)) domain
  and labels2 l = 
    Array.map (fun (x,y) -> Printf.sprintf "%s %s/%s" l (B.Chipseq.string_of_sample x) (B.Chipseq.string_of_sample y)) domain2
  in
  Array.enum annot
  /@ Location.(fun a -> Array.concat [
    [| a.loc.chr ; string_of_int a.loc.st ; string_of_int a.loc.ed |] ;
    Array.map string_of_int a.counts ;
    Array.map string_of_float a.coverage ;
    Array.map string_of_float a.pvalue 
  ])
  |> Enum.(append (singleton (
       Array.concat [
	 [| "chrom" ; "start" ; "end" |] ;
	 labels "counts" ;
	 labels "coverage" ;
	 labels2 "pvalue" 
       ])))
  |> Enum.map Tsv.string_of_line
  |> File.write_lines output_path

let i (_,samples,_) x = Array.findi (( = ) x) samples
let i2 (_,_,design) x y = Array.findi (( = ) (x,y)) design




















