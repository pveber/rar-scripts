open Batteries

let uniform ?seed ~chrom_size loc items n = 
  let rng = Gsl_rng.make (Gsl_rng.default ()) in
  Option.may (fun seed -> Gsl_rng.set rng (Nativeint.of_int seed)) seed ;
  let pick_chr = 
    let p = Array.map (fun x -> float (snd x)) chrom_size in
    let d = Gsl_randist.discrete_preproc p in
    fun () -> chrom_size.(Gsl_randist.discrete rng d)
  in
  (1 -- n) 
  /@ (fun _ -> 
    let loc = loc items.(Gsl_rng.uniform_int rng (Array.length items)) in
    let width = Location.size loc in
    let chr, chr_length = pick_chr () in 
    let k = Gsl_rng.uniform_int rng (chr_length - width) in 
    Location.make chr k (k + width - 1))





















