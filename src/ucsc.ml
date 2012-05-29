open Batteries
open Printf

type genome = [ `mm9 | `hg18 ]
let string_of_genome = function
    `mm9 -> "mm9"
  | `hg18 -> "hg18"

let chrom_size_cmd1 org = sprintf "\
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
-e 'select chrom,size from chromInfo;' %s" (string_of_genome org)
let chrom_size_cmd2 path = sprintf "\
gawk -F'\t' '{printf \"%%s\\t%%s\\n\", $1,$2}' >> %s" path

let fetch_chrom_size org path =
  Sle.sh "rm -f %s" path ;
  Sle.pipefail (chrom_size_cmd1 org) (chrom_size_cmd2 path)

let chrom_size_of_file path =
  Tsv.enum path
  /@ (fun l -> l.(0), int_of_string l.(1))





















