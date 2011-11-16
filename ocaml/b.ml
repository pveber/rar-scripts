(* #use "topfind" *)
(* #require "extlib" *)

open Printf
open Genome
open Oregon
open Sle.Infix

let ncores = Array.length (Sle.shout "cat /proc/cpuinfo | grep '^processor'")

(** available memory in MB *)
let mem = int_of_string (Sle.SleString.split_noeps ' ' (Sle.shout "free -mt").(1)).(1);;
let preview_mode = mem < 4000

module Setting = struct
  let ncores = ncores
  let org = `mouse

  type condition = [
    `t0
  | `dmso of int
  | `atra of int 
  ]
  and condition_type = condition

  let string_of_condition = function 
      `t0 -> "t0"
    | `dmso t -> "DMSO-t" ^ (string_of_int t)
    | `atra t -> "ATRA-t" ^ (string_of_int t)

  module RnaseqData = struct
    type sample = condition * [ `polyA | `riboM ]
    type condition = condition_type

    let string_of_sample (condition,technique) = 
      (string_of_condition condition) ^ "-" ^ (match technique with `polyA -> "polyA" | `riboM -> "rbm")

    let samples = [ (*(`t0,      `riboM) ; 
		    (`dmso 24, `riboM) ; 
		    (`dmso 48, `riboM) ; 
		    (`atra 12, `riboM) ; 
		    (`atra 24, `riboM) ; 
		    (`atra 36, `riboM) ; 
		    (`atra 48, `riboM) ; *)

		    (`t0,      `polyA) ;
		    (`atra 6,  `polyA) ;
		    (`atra 12, `polyA) ;
		    (`atra 24, `polyA) ;
		    (`atra 36, `polyA) ;
		    (`atra 48, `polyA) ;
		    (`dmso 24, `polyA) ;
		    (`dmso 48, `polyA)
    ]

    let read_length _ = 72

    let replicates = function x -> [ x ]

    type channel = sample

    let channels = replicates 

    let path = function
      | `t0,      `riboM  -> "data-gb/rnaseq-20100616/s_1_sequence"
      | `dmso 24, `riboM -> "data-gb/rnaseq-20100616/s_2_sequence"
      | `dmso 48, `riboM -> "data-gb/rnaseq-20100616/s_3_sequence"
      | `atra 12, `riboM -> "data-gb/rnaseq-20100616/s_4_sequence"
      | `atra 24, `riboM -> "data-gb/rnaseq-20100616/s_5_sequence"
      | `atra 36, `riboM -> "data-gb/rnaseq-20100616/s_6_sequence"
      | `atra 48, `riboM -> "data-gb/rnaseq-20100616/s_7_sequence"
      | `t0,      `polyA -> "data-gb/rnaseq-20110211/s_1_sequence.txt"
      | `dmso 24, `polyA -> "data-gb/rnaseq-20110211/s_3_sequence.txt"
      | `dmso 48, `polyA -> "data-gb/rnaseq-20110211/s_5_sequence.txt"
      | `atra  6, `polyA -> "data-gb/rnaseq-20110211/s_2_sequence.txt"
      | `atra 12, `polyA -> "data-gb/rnaseq-20110211/s_7_sequence.txt"
      | `atra 24, `polyA -> "data-gb/rnaseq-20110211/s_4_sequence.txt"
      | `atra 36, `polyA -> "data-gb/rnaseq-20110211/s_8_sequence.txt"
      | `atra 48, `polyA -> "data-gb/rnaseq-20110211/s_6_sequence.txt"
      | _ -> assert false

    let fastq_type _ = `fastq `solexa

    let ribosomal_transcripts_fasta = [|
      { Fasta.comment = "NR_030686.1" ;
	Fasta.seq = "gtctacggccataccaccctgaacgcgcccgatctcgtctgatctcggaagctaagcagggtcgggcctggttagtacttggatgggagaccgcctgggaataccgggtgctgtaggcttt" } ;
      { Fasta.comment = "NR_003280.1" ;
	Fasta.seq = "cgactcttagcggtggatcactcggctcgtgcgtcgatgaagaacgcagcgctagctgcgagaattaatgtgaattgcaggacacattgatcatcgacacttcgaacgcacttgcggccccgggttcctcccggggctacgcctgtctgagcgtcgct" } ;
      { Fasta.comment = "NR_003278.1" ;
	Fasta.seq = "acctggttgatcctgccaggtagcatatgcttgtctcaaagattaagccatgcatgtctaagtacgcacggccggtacagtgaaactgcgaatggctcattaaatcagttatggttcctttggtcgctcgctcctctcctacttggataactgtggtaattctagagctaatacatgccgacgggcgctgaccccccttcccggggggggatgcgtgcatttatcagatcaaaaccaacccggtgagctccctcccggctccggccgggggtcgggcgccggcggcttggtgactctagataacctcgggccgatcgcacgccccccgtggcggcgacgacccattcgaacgtctgccctatcaactttcgatggtagtcgccgtgcctaccatggtgaccacgggtgacggggaatcagggttcgattccggagagggagcctgagaaacggctaccacatccaaggaaggcagcaggcgcgcaaattacccactcccgacccggggaggtagtgacgaaaaataacaatacaggactctttcgaggccctgtaattggaatgagtccactttaaatcctttaacgaggatccattggagggcaagtctggtgccagcagccgcggtaattccagctccaatagcgtatattaaagttgctgcagttaaaaagctcgtagttggatcttgggagcgggcgggcggtccgccgcgaggcgagtcaccgcccgtccccgccccttgcctctcggcgccccctcgatgctcttagctgagtgtcccgcggggcccgaagcgtttactttgaaaaaattagagtgttcaaagcaggcccgagccgcctggataccgcagctaggaataatggaataggaccgcggttctattttgttggttttcggaactgaggccatgattaagagggacggccgggggcattcgtattgcgccgctagaggtgaaattcttggaccggcgcaagacggaccagagcgaaagcatttgccaagaatgttttcattaatcaagaacgaaagtcggaggttcgaagacgatcagataccgtcgtagttccgaccataaacgatgccgactggcgatgcggcggcgttattcccatgacccgccgggcagcttccgggaaaccaaagtctttgggttccggggggagtatggttgcaaagctgaaacttaaaggaattgacggaagggcaccaccaggagtgggcctgcggcttaatttgactcaacacgggaaacctcacccggcccggacacggacaggattgacagattgatagctctttctcgattccgtgggtggtggtgcatggccgttcttagttggtggagcgatttgtctggttaattccgataacgaacgagactctggcatgctaactagttacgcgacccccgagcggtcggcgtcccccaacttcttagagggacaagtggcgttcagccacccgagattgagcaataacaggtctgtgatgcccttagatgtccggggctgcacgcgcgctacactgactggctcagcgtgtgcctaccctgcgccggcaggcgcgggtaacccgttgaaccccattcgtgatggggatcggggattgcaattattccccatgaacgaggaattcccagtaagtgcgggtcataagcttgcgttgattaagtccctgccctttgtacacaccgcccgtcgctactaccgattggatggtttagtgaggccctcggatcggccccgccggggtcggcccacggccctggcggagcgctgagaagacggtcgaacttgactatctagaggaagtaaaagtcgtaacaaggtttccgtaggtgaacctgcggaaggatcattaa" } ;
      { Fasta.comment = "NR_003279.1" ;
	Fasta.seq = "cgcgacctcagatcagacgtggcgacccgctgaatttaagcatattagtcagcggaggaaaagaaactaaccaggattccctcagtaacggcgagtgaacagggaagagcccagcgccgaatccccgccgcgcgtcgcggcgtgggaaatgtggcgtacggaagacccactccccggcgccgctcgtggggggcccaagtccttctgatcgaggcccagcccgtggacggtgtgaggccggtagcggccccggcgcgccgggctcgggtcttcccggagtcgggttgcttgggaatgcagcccaaagcgggtggtaaactccatctaaggctaaataccggcacgagaccgatagtcaacaagtaccgtaagggaaagttgaaaagaactttgaagagagagttcaagagggcgtgaaaccgttaagaggtaaacgggtggggtccgcgcagtccgcccggaggattcaacccggcggcgcgcgtccggccgtgcccggtggtcccggcggatctttcccgctccccgttcctcccgacccctccacccgcgcgtcgttcccctcttcctccccgcgtccggcgcctccggcggcgggcgcggggggtggtgtggtggtggcgcgcgggcggggccgggggtggggtcggcgggggaccgcccccggccggcgaccggccgccgccgggcgcacttccaccgtggcggtgcgccgcgaccggctccgggacggccgggaaggcccggtggggaaggtggctcggggggggcggcgcgtctcagggcgcgccgaaccacctcaccccgagtgttacagccctccggccgcgctttcgccgaatcccggggccgaggaagccagatacccgtcgccgcgctctccctctccccccgtccgcctcccgggcgggcgtgggggtgggggccgggccgcccctcccacggcgcgaccgctctcccacccccctccgtcgcctctctcggggcccggtggggggcggggcggactgtccccagtgcgccccgggcgtcgtcgcgccgtcgggtcccggggggaccgtcggtcacgcgtctcccgacgaagccgagcgcacggggtcggcggcgatgtcggctacccacccgacccgtcttgaaacacggaccaaggagtctaacgcgtgcgcgagtcaggggctcgtccgaaagccgccgtggcgcaatgaaggtgaagggccccgcccgggggcccgaggtgggatcccgaggcctctccagtccgccgagggcgcaccaccggcccgtctcgcccgccgcgccggggaggtggagcacgagcgtacgcgttaggacccgaaagatggtgaactatgcttgggcagggcgaagccagaggaaactctggtggaggtccgtagcggtcctgacgtgcaaatcggtcgtccgacctgggtataggggcgaaagactaatcgaaccatctagtagctggttccctccgaagtttccctcaggatagctggcgctctcgctcccgacgtacgcagttttatccggtaaagcgaatgattagaggtcttggggccgaaacgatctcaacctattctcaaactttaaatgggtaagaagcccggctcgctggcgtggagccgggcgtggaatgcgagtgcctagtgggccacttttggtaagcagaactggcgctgcgggatgaaccgaacgccgggttaaggcgcccgatgccgacgctcatcagaccccagaaaaggtgttggttgatatagacagcaggacggtggccatggaagtcggaatccgctaaggagtgtgtaacaactcacctgccgaatcaactagccctgaaaatggatggcgctggagcgtcgggcccatacccggccgtcgccgcagtcggaacggaacgggacgggagcggccgcgggtgcgcgtctctcggggtcgggggtgcgtggcgggggcccgtcccccgcctcccctccgcgcgccgggttcgcccccgcggcgtcgggccccgcggagcctacgccgcgacgagtaggagggccgctgcggtgagccttgaagcctagggcgcgggcccgggtggagccgccgcaggtgcagatcttggtggtagtagcaaatattcaaacgagaactttgaaggccgaagtggagaagggttccatgtgaacagcagttgaacatgggtcagtcggtcctgagagatgggcgagtgccgttccgaagggacgggcgatggcctccgttgccctcggccgatcgaaagggagtcgggttcagatccccgaatccggagtggcggagatgggcgccgcgaggccagtgcggtaacgcgaccgatcccggagaagccggcgggaggcctcggggagagttctcttttctttgtgaagggcagggcgccctggaatgggttcgccccgagagaggggcccgtgccttggaaagcgtcgcggttccggcggcgtccggtgagctctcgctggcccttgaaaatccgggggagagggtgtaaatctcgcgccgggccgtacccatatccgcagcaggtctccaaggtgaacagcctctggcatgttggaacaatgtaggtaagggaagtcggcaagccggatccgtaacttcgggataaggattggctctaagggctgggtcggtcgggctggggcgcgaagcggggctgggcgcgcgccgcggctggacgaggcgccgccgccctctcccacgtccggggagaccccccgtcctttccgcccgggcccgccctcccctcttccccgcggggccccgtcgtcccccgcgtcgtcgccacctctcttcccccctccttcttcccgtcggggggcgggtcgggggtcggcgcgcggcgcgggctccggggcggcgggtccaaccccgcgggggttccggagcgggaggaaccagcggtccccggtggggcggggggcccggacactcggggggccggcggcggcggcgactctggacgcgagccgggcccttcccgtggatcgcctcagctgcggcgggcgtcgcggccgctcccggggagcccggcgggtgccggcgcgggtcccctccccgcggggcctcgctccacccccccatcgcctctcccgaggtgcgtggcgggggcgggcgggcgtgtcccgcgcgtgtggggggaacctccgcgtcggtgttcccccgccgggtccgccccccgggccgcggttttccgcgcggcgcccccgcctcggccggcgcctagcagccgacttagaactggtgcggaccaggggaatccgactgtttaattaaaacaaagcatcgcgaaggcccgcggcgggtgttgacgcgatgtgatttctgcccagtgctctgaatgtcaaagtgaagaaattcaatgaagcgcgggtaaacggcgggagtaactatgactctcttaaggtagccaaatgcctcgtcatctaattagtgacgcgcatgaatggatgaacgagattcccactgtccctacctactatccagcgaaaccacagccaagggaacgggcttggcggaatcagcggggaaagaagaccctgttgagcttgactctagtctggcacggtgaagagacatgagaggtgtagaataagtgggaggcccccggcgcccggccccgtcctcgcgtcggggtcggggcacgccggcctcgcgggccgccggtgaaataccactactctcatcgttttttcactgacccggtgaggcgggggggcgagccccgaggggctctcgcttctggcgccaagcgtccgtcccgcgcgtgcgggcgggcgcgacccgctccggggacagtgccaggtggggagtttgactggggcggtacacctgtcaaacggtaacgcaggtgtcctaaggcgagctcagggaggacagaaacctcccgtggagcagaagggcaaaagctcgcttgatcttgattttcagtacgaatacagaccgtgaaagcggggcctcacgatccttctgaccttttgggttttaagcaggaggtgtcagaaaagttaccacagggataactggcttgtggcggccaagcgttcatagcgacgtcgctttttgatccttcgatgtcggctcttcctatcattgtgaagcagaattcaccaagcgttggattgttcacccactaatagggaacgtgagctgggtttagaccgtcgtgagacaggttagttttaccctactgatgatgtgttgttgccatggtaatcctgctcagtacgagaggaaccgcaggttcagacatttggtgtatgtgcttggctgaggagccaatggggcgaagctaccatctgtgggattatgactgaacgcctctaagtcagaatcccgcccaggcggaacgatacggcagcgccgaaggagcctcggttggccccggatagccgggtccccgtccgtccccgctcggcggggtccccgcgtcgtccccgcggcggcgcggggtctccccccgccgggcgtcgggaccggggtccggtgcggagagccgttcgtcttgggaaacggggtgcggccggaaagggggccgccctctcgcccgtcacgttgaacgcacgttcgtgtggaacctggcgctaaaccattcgtagacgacctgcttctgggtcggggtttcgtacgtagcagagcagctccctcgctgcgatctattgaaagtcagccctcgacacaagggtttgtctctgcgggctttc" } 
    |]

    let ribosomal_transcripts = Target.F.make 
      (object
	 method id = "Gerard.ribosomal_transcripts"
	 method deps = []
	 method build path =
	   Fasta.to_file path ribosomal_transcripts_fasta
	 method ty = `fasta
       end)

    let pretreatment = function 
      | (_, `riboM) -> fun x ->
	Ngs_utils.filter_reads 
	  x 
	  (Bowtie.index 
	     ~packed:true
	     ribosomal_transcripts)
      | (_, `polyA) -> fun x -> x

    let comp x y = 
      (string_of_condition x) ^ "/" ^ (string_of_condition y),
      [ (x, `polyA) ],
      [ (y, `polyA) ]

    let comp_t0 x = 
      (string_of_condition `t0) ^ "/" ^ (string_of_condition x),
      [ (`t0, `polyA) ; (`dmso 24, `polyA) ; (`dmso 48, `polyA) ],
      [ (x, `polyA) ]

    let design = [
      comp_t0 (`atra 6);
      comp_t0 (`atra 12);
      comp_t0 (`atra 24);
      comp_t0 (`atra 36);
      comp_t0 (`atra 48);
      comp_t0 (`dmso 24);
      comp_t0 (`dmso 48);
      comp (`atra 6) (`atra 12);
      comp (`atra 12) (`atra 24);
      comp (`atra 24) (`atra 36);
      comp (`atra 36) (`atra 48);
    ]	
  end

  module ChipseqData = struct
    type sample = [
	`F9_ATRA_panRAR_1
      | `F9_ATRA_panRAR_2
      | `Input_1
      | `F9_WT_panRAR_1
      | `F9_WT_panRXR_1
      | `F9_ATRA24_panRAR_1
      | `F9_ATRA24_panRXR_1
      | `F9_ATRA48_panRAR_1
      | `F9_ATRA48_panRXR_1
      | `F9_ATRA2_panRXR_1
      | `F9_WT_PolII_1
      | `F9_ATRA2_PolII_1
      | `F9_ATRA24_PolII_1
      | `F9_ATRA48_PolII_1
      | `Mendoza_input
      | `Mendoza_ATRA2_RARg
      | `Mendoza_ATRA6_RARg
      | `Mendoza_ATRA24_RARg
      | `Mendoza_ATRA48_RARg
      | `Mendoza_EtOH48_RARg
      | `Mendoza_ATRA2_RXRa
      | `Mendoza_ATRA6_RXRa
      | `Mendoza_ATRA24_RXRa
      | `Mendoza_ATRA48_RXRa
      | `Mendoza_EtOH48_RXRa
    ]
    type condition = condition_type

    let string_of_sample = function
	`F9_ATRA_panRAR_1 -> "F9-ATRA-panRAR-1"
      | `F9_ATRA_panRAR_2 -> "F9-ATRA-panRAR-2"
      | `Input_1 -> "input-1"
      | `F9_WT_panRAR_1 -> "F9-WT-panRAR-1"
      | `F9_WT_panRXR_1 -> "F9-WT-panRXR-1"
      | `F9_ATRA24_panRAR_1 -> "F9-ATRA24-panRAR-1"
      | `F9_ATRA24_panRXR_1 -> "F9-ATRA24-panRXR-1"
      | `F9_ATRA48_panRAR_1 -> "F9-ATRA48-panRAR-1"
      | `F9_ATRA48_panRXR_1 -> "F9-ATRA48-panRXR-1"
      | `F9_ATRA2_panRXR_1 -> "F9-ATRA2-panRXR-1"
      | `F9_WT_PolII_1 -> "F9-WT-PolII-1"
      | `F9_ATRA2_PolII_1 -> "F9-ATRA2-PolII-1"
      | `F9_ATRA24_PolII_1 -> "F9-ATRA24-PolII-1"
      | `F9_ATRA48_PolII_1 -> "F9-ATRA48-PolII-1"
      | `Mendoza_input -> "Mendoza-WCE"
      | `Mendoza_ATRA2_RARg -> "Mendoza-ATRA2-RARg"
      | `Mendoza_ATRA6_RARg -> "Mendoza-ATRA6-RARg"
      | `Mendoza_ATRA24_RARg -> "Mendoza-ATRA24-RARg"
      | `Mendoza_ATRA48_RARg -> "Mendoza-ATRA48-RARg"
      | `Mendoza_EtOH48_RARg -> "Mendoza-EtOH48-RARg"
      | `Mendoza_ATRA2_RXRa -> "Mendoza-ATRA2-RXRa"
      | `Mendoza_ATRA6_RXRa -> "Mendoza-ATRA6-RXRa"
      | `Mendoza_ATRA24_RXRa -> "Mendoza-ATRA24-RXRa"
      | `Mendoza_ATRA48_RXRa -> "Mendoza-ATRA48-RXRa"
      | `Mendoza_EtOH48_RXRa -> "Mendoza-EtOH48-RXRa"


    let samples = [ `F9_ATRA_panRAR_1 ; `F9_ATRA_panRAR_2 ; `Input_1 ;       
		    `F9_WT_panRAR_1 ; `F9_WT_panRXR_1 ; `F9_ATRA24_panRAR_1 ;
		    `F9_ATRA24_panRXR_1 ; `F9_ATRA48_panRAR_1 ; `F9_ATRA48_panRXR_1 ;
		    `F9_ATRA2_panRXR_1 ; `F9_WT_PolII_1 ; `F9_ATRA2_PolII_1 ; 
		    `F9_ATRA24_PolII_1 ; `F9_ATRA48_PolII_1 ;
		    `Mendoza_input ; `Mendoza_ATRA2_RARg ; `Mendoza_ATRA6_RARg ; 
		    `Mendoza_ATRA24_RARg ; `Mendoza_ATRA48_RARg ; `Mendoza_EtOH48_RARg ;
		    `Mendoza_ATRA2_RXRa ; `Mendoza_ATRA6_RXRa ; 
		    `Mendoza_ATRA24_RXRa ; `Mendoza_ATRA48_RXRa ; `Mendoza_EtOH48_RXRa ]

    type channel = sample

    let channels x = [ x ]

    let path = function
	`F9_ATRA_panRAR_1 -> "data-gb/chipseq-20100709/s_1_sequence"
      | `F9_ATRA_panRAR_2 -> "data-gb/chipseq-20100709/s_2_sequence"
      | `Input_1 -> "data-gb/chipseq-20100709/s_3_sequence"
      | `F9_WT_panRAR_1 -> "data-gb/chipseq-20101202/s_1_sequence.txt"
      | `F9_WT_panRXR_1 -> "data-gb/chipseq-20101202/s_2_sequence.txt"
      | `F9_ATRA24_panRAR_1 -> "data-gb/chipseq-20101202/s_3_sequence.txt"
      | `F9_ATRA24_panRXR_1 -> "data-gb/chipseq-20101202/s_4_sequence.txt"
      | `F9_ATRA48_panRAR_1 -> "data-gb/chipseq-20101202/s_5_sequence.txt"
      | `F9_ATRA48_panRXR_1 -> "data-gb/chipseq-20101202/s_6_sequence.txt"
      | `F9_ATRA2_panRXR_1 -> "data-gb/chipseq-20101202/s_7_sequence.txt"
      | `F9_WT_PolII_1 -> "data-gb/chipseq-20110413/s_3_sequence.txt"
      | `F9_ATRA2_PolII_1 -> "data-gb/chipseq-20110413/s_4_sequence.txt"
      | `F9_ATRA24_PolII_1 -> "data-gb/chipseq-20110413/s_5_sequence.txt"
      | `F9_ATRA48_PolII_1 -> "data-gb/chipseq-20110413/s_6_sequence.txt"
      | `Mendoza_input -> "data-mendoza/input.fastq"
      | `Mendoza_ATRA2_RARg -> "data-mendoza/ATRA2_RARg.fastq"
      | `Mendoza_ATRA6_RARg -> "data-mendoza/ATRA6_RARg.fastq"
      | `Mendoza_ATRA24_RARg -> "data-mendoza/ATRA24_RARg.fastq"
      | `Mendoza_ATRA48_RARg -> "data-mendoza/ATRA48_RARg.fastq"
      | `Mendoza_EtOH48_RARg -> "data-mendoza/EtOH48_RARg.fastq"
      | `Mendoza_ATRA2_RXRa -> "data-mendoza/ATRA2_RXRa.fastq"
      | `Mendoza_ATRA6_RXRa -> "data-mendoza/ATRA6_RXRa.fastq"
      | `Mendoza_ATRA24_RXRa -> "data-mendoza/ATRA24_RXRa.fastq"
      | `Mendoza_ATRA48_RXRa -> "data-mendoza/ATRA48_RXRa.fastq"
      | `Mendoza_EtOH48_RXRa -> "data-mendoza/EtOH48_RXRa.fastq"


    let fastq_type = function
      | `Mendoza_input -> `fastq `phred33
      | `Mendoza_ATRA2_RARg -> `fastq `phred33
      | `Mendoza_ATRA6_RARg -> `fastq `phred33
      | `Mendoza_ATRA24_RARg -> `fastq `phred33
      | `Mendoza_ATRA48_RARg -> `fastq `phred33
      | `Mendoza_EtOH48_RARg -> `fastq `phred33
      | `Mendoza_ATRA2_RXRa -> `fastq `phred33
      | `Mendoza_ATRA6_RXRa -> `fastq `phred33
      | `Mendoza_ATRA24_RXRa -> `fastq `phred33
      | `Mendoza_ATRA48_RXRa -> `fastq `phred33
      | `Mendoza_EtOH48_RXRa -> `fastq `phred33
      | _ -> `fastq `solexa

    let design = [ `F9_ATRA_panRAR_1, `Input_1 ; `F9_ATRA_panRAR_2, `Input_1 ;
		   `F9_WT_panRAR_1, `Input_1 ; `F9_WT_panRXR_1, `Input_1 ; 
		   `F9_ATRA24_panRAR_1, `Input_1 ; `F9_ATRA24_panRXR_1, `Input_1 ;
		   `F9_ATRA48_panRAR_1, `Input_1 ; `F9_ATRA48_panRXR_1, `Input_1 ;
		   `F9_ATRA2_panRXR_1, `Input_1 ; `F9_WT_PolII_1, `Input_1 ; 
		   `F9_ATRA2_PolII_1, `Input_1 ; `F9_ATRA24_PolII_1, `Input_1 ; 
		   `F9_ATRA48_PolII_1, `Input_1 ;
		   `Mendoza_ATRA2_RARg, `Mendoza_input ; `Mendoza_ATRA6_RARg, `Mendoza_input ;
		   `Mendoza_ATRA24_RARg, `Mendoza_input ; `Mendoza_ATRA48_RARg, `Mendoza_input ;
		   `Mendoza_EtOH48_RARg, `Mendoza_input ;
		   `Mendoza_ATRA2_RXRa, `Mendoza_input ; `Mendoza_ATRA6_RXRa, `Mendoza_input ;
		   `Mendoza_ATRA24_RXRa, `Mendoza_input ; `Mendoza_ATRA48_RXRa, `Mendoza_input ;
		   `Mendoza_EtOH48_RXRa, `Mendoza_input ]

    let pretreatment y x = x

    let macs_mfold_default = 16
    let macs_mfold = [ `Input_1, 10 ]
  end 

  module Motif_library = struct
    open Target.Infix
    open Motif_library

    let balmer_hexamer = `sequence [
      `base (1.04, -1.8,-0.01,-1.78) ;
      `base (-2.47,-3.34,1.27,-1.87) ;
      `base (-1.7,-3.05,0.94,0.16) ;
      `base (-2.88,-1.3,-0.69,1.22) ;
      `base (-3.17,1.09,-0.85,-1.38) ;
      `base (1.35,-3.34,-1.36,-2.07) ;
    ]

    let balmer_dr125 = 
      PSSM (`sequence [ 
	      balmer_hexamer ; 
	      `disjunction [ 
		`gap (5,5) ; 
		`gap (2,2) ; 
		`gap (1,1) ] ; 
	      balmer_hexamer ],
	    "DR1-2-5")

    let balmer_dr0125 = 
      PSSM (`sequence [ 
	      balmer_hexamer ; 
	      `disjunction [ 
		`gap (0,0) ; 
		`gap (5,5) ; 
		`gap (2,2) ; 
		`gap (1,1) ] ; 
	      balmer_hexamer ],
	    "DR0-1-2-5")

    let balmer_dr0_5 = `sequence [
      balmer_hexamer ;
      `gap (0,5) ; 
      balmer_hexamer
    ]

    let balmer_drn i =
      PSSM (`sequence [ balmer_hexamer ; `gap (i,i) ; balmer_hexamer ],
	    sprintf "DR%d" i)
	
    let balmer_dr0_n i =
      PSSM (`sequence [ balmer_hexamer ; `gap (0,i) ; balmer_hexamer ],
	    sprintf "BalmerDR0-%d" i)

    let balmer_dtr0_n i =
      PSSM (`sequence [ balmer_hexamer ; `gap (0,i) ; balmer_hexamer ; `gap (0,i) ; balmer_hexamer ],
	    sprintf "BalmerDTR0-%d" i)

    let composite_trn i =
      PSSM (`disjunction [
	      `sequence [ balmer_hexamer ; `gap (i,i) ; balmer_hexamer ; `disjunction [ `gap (1,1) ;`gap (2,2) ;`gap (5,5) ] ; balmer_hexamer ] ;
	      `sequence [ balmer_hexamer ; `disjunction [ `gap (1,1) ;`gap (2,2) ;`gap (5,5) ] ; balmer_hexamer ; `gap (i,i) ; balmer_hexamer ]
	    ],
	    sprintf "BalmerCompositeTR-%d" i)

    let balmer_dqr0_n i =
      PSSM (`sequence [ balmer_hexamer ; `gap (0,i) ; balmer_hexamer ; `gap (0,i) ; balmer_hexamer ; `gap (0,i) ; balmer_hexamer ],
	    sprintf "BalmerDQR0-%d" i)

    let balmer_irn i =
      PSSM (`sequence [ balmer_hexamer ; `gap (i,i) ; Motif.PSSM.reverse_complement balmer_hexamer ],
	    sprintf "IR%d" i)

    let balmer_ern i =
      PSSM (`sequence [ Motif.PSSM.reverse_complement balmer_hexamer ; `gap (i,i) ; balmer_hexamer ],
	    sprintf "ER%d" i)
	
    let value = Array.concat [
      Array.init 6 balmer_drn ;
      Array.init 6 balmer_ern ;
      Array.init 6 balmer_irn ;
      [| balmer_dr125 ; balmer_dr0125 |]
    ]
    let deps  = []
  end
end

module B = Bench.Make(Setting);;

let rep = if preview_mode then "benoit_preview" else "benoit"
let _ = GzmConfig.load ".guizminrc.benoit" 
let login = GzmConfig.export_login ()
let passwd = GzmConfig.export_passwd ()
let url = sprintf "http://%s:%s@cgmc.univ-lyon1.fr/FGNRseq/%s" login passwd rep
let host = GzmConfig.host_ssh_url () ^ "/" ^ rep
let website = B.Output.website ~url ~host;;

open Motif_library
open Setting.Motif_library

let jaspar_tandems motif = 
  motif_tandems ~mingap:0 ~maxgap:50 motif (Motif_library.of_jaspar_collection Jaspar.core)


include B



let motif_study peaks =
  let peakseq = Ucsc.fasta_of_bed Genome.id peaks in 
  let ctrlseq = Fasta.enumerate (Ucsc.control_sequences Genome.id peakseq) in
  (object
     method roc_balmer125 = 
       roc_curve balmer_dr125 peakseq ctrlseq

     method roc_balmer0_5 = 
       roc_curve 
	 (PSSM (balmer_dr0_5,"BalmerDR0-5"))
	 peakseq ctrlseq

     method roc_bundle_drn = 
       roc_bundle 
	 (Array.init 6 balmer_drn)
	 peakseq ctrlseq

     method drn_ranking = 
       rank
	 (Array.init 15 balmer_drn)
	 peakseq ctrlseq

     method irn_ranking = 
       rank
	 (Array.init 15 balmer_irn)
	 peakseq ctrlseq

     method ern_ranking = 
       rank
	 (Array.init 15 balmer_ern)
	 peakseq ctrlseq

     method dtr0_n_ranking = 
       rank
	 (Array.init 15 balmer_dtr0_n)
	 peakseq ctrlseq

     method composite_trn_ranking = 
       rank
	 (Array.init 15 composite_trn)
	 peakseq ctrlseq

     method dqr0_n_ranking = 
       rank
	 (Array.init 15 balmer_dqr0_n)
	 peakseq ctrlseq

     method dr0_n_ranking = 
       rank
	 (Array.init 15 balmer_dr0_n)
	 peakseq ctrlseq

     method dr0_10_roc = 
       roc_curve
	 (balmer_dr0_n 10)
	 peakseq ctrlseq

     method jaspar_ranking = 
       rank
	 (Motif_library.of_jaspar_collection Jaspar.core)
	 peakseq ctrlseq
	 
     method jaspar_ranking_shuffle_test = 
       rank
	 (Motif_library.of_jaspar_collection Jaspar.core)
	 peakseq (Fasta.shuffle peakseq)

     method jaspar_tandems_ranking motif =
       rank (jaspar_tandems motif) peakseq ctrlseq

     method meme_hexamer =
       let prior = MemeSuite.psp_gen ~revcomp:true peakseq ctrlseq in
       MemeSuite.meme ~minw:6 ~maxw:6 ~nmotifs:5 ~revcomp:true ~mode:`zoops ~psp:prior peakseq
     method meme =
       let prior = MemeSuite.psp_gen ~revcomp:true peakseq ctrlseq in
       MemeSuite.meme ~minw:6 ~maxw:20 ~nmotifs:10 ~revcomp:true ~mode:`zoops ~psp:prior peakseq
   end)


let peaks = Chipseq.best_macs_peaks 10000 `F9_ATRA_panRAR_1
let peaksummits = Macs.locations_around_summit ~radius:100 peaks

let ms = motif_study peaks
let mss = motif_study peaksummits


module M(P : sig end) = struct
  module L = Labs.Make(B)
  open L
  open Batteries

  let validated_chipseq_samples = 
    [ `F9_ATRA_panRAR_1 ;
      `F9_WT_panRAR_1 ; `F9_WT_panRXR_1 ; `F9_ATRA24_panRAR_1 ;
      `F9_ATRA24_panRXR_1 ; `F9_ATRA48_panRAR_1 ; `F9_ATRA48_panRXR_1 ;
      `F9_ATRA2_panRXR_1 ]

  let regions_bed = binding_loci_tsv validated_chipseq_samples
  let read_counts =
    Sle.hfun_make
      (fun x -> (Ngs_utils.read_counts_of_location_file (Chipseq.bowtie_wodup x) regions_bed)#value)
      Chipseq.samples
  let library_size = Sle.hfun_make (fun x -> (Sam.nbmappings (Chipseq.bowtie_wodup x))#value) Chipseq.samples

  let pmt_control ~treatment ~control ~treatment_size ~control_size = 
    Pmt.log_poisson_margin control treatment control_size treatment_size


  let peak_islands = (binding_loci validated_chipseq_samples)#value |> Array.of_list

  let regions = Array.map (fun l -> let x,_,_ = peak_average l in Location.(make x.chr (x.st - 200) (x.ed + 200))) peak_islands

  let profiles_against_input = 
    Array.mapi
      (fun i _ -> 
	 let control = (`Input_1 --> read_counts).(i)
	 and control_size = `Input_1 --> library_size in
	 Sle.hfun_make 
	   (fun x -> -. pmt_control ~treatment:(x --> read_counts).(i) ~control ~treatment_size:(x --> library_size) ~control_size)
	   validated_chipseq_samples)
      regions

  let intensity = 
    Array.mapi
      (fun i _ -> 
	 Sle.hfun_make 
	   (fun x -> float (x --> read_counts).(i) /. float (x --> library_size) *. 1e6)
	   validated_chipseq_samples)
      regions

  let binding_loci_pmt_tsv path = 
    let samples = Array.of_list validated_chipseq_samples in
    (0 --^ Array.length regions) 
    |> Enum.map (fun i -> Location.(let l = regions.(i) in
				    Array.concat [
				      [| l.chr ; string_of_int l.st ; string_of_int l.ed |] ;
				      Array.map (fun x -> sprintf "%g" (x --> profiles_against_input.(i))) samples ;
				      Array.map (fun x -> sprintf "%g" (x --> intensity.(i))) samples ]))
    |> Enum.append (Enum.singleton (
		      Array.concat [
			[| "chr" ; "start" ; "end" |] ;
			Array.map (fun x -> (Chipseq.string_of_sample x) ^ " pval") samples ;
			Array.map (fun x -> (Chipseq.string_of_sample x) ^ " intens") samples ] ))
    |> Enum.map Tsv.string_of_line
    |> File.write_lines path

  let map f = Array.mapi (fun i _ -> f i) regions

  let fold mask f init = Array.fold_lefti (fun accu i b -> if b then f i accu else accu) init mask
    
  let nbregions mask = fold mask (fun _ -> (+) 1) 0

  let filter mask t = Sle.SleArray.filteri (fun i _ -> mask.(i)) t

  let _and m1 m2 = Array.mapi (fun i _ -> m1.(i) && m2.(i)) m1

  let rar_rxr_pairs = [
    `F9_WT_panRAR_1, `F9_WT_panRXR_1 ;
    `F9_ATRA_panRAR_1, `F9_ATRA2_panRXR_1 ;
    `F9_ATRA24_panRAR_1, `F9_ATRA24_panRXR_1 ;
    `F9_ATRA48_panRAR_1, `F9_ATRA48_panRXR_1 ;
  ]

  let binding theta x i = 
    (x --> profiles_against_input.(i)) > theta

  let binding_in_some_condition theta i = 
    List.exists (fun x -> binding theta x i) validated_chipseq_samples

  let rar_binding_in_some_condition theta i = 
    List.exists (fun (rar,rxr) -> binding theta rar i) rar_rxr_pairs

  let rxr_binding_in_some_condition theta i = 
    List.exists (fun (rar,rxr) -> binding theta rxr i) rar_rxr_pairs

  let confirmed_binding_in_some_condition theta i = 
    List.exists (fun (rar,rxr) -> binding theta rar i && binding theta rxr i) rar_rxr_pairs

  let ( -| ) x y = Array.mapi (fun i _ -> x.(i) && (not y.(i))) x

  let etages f = 
    Array.init 10 
      (fun i -> 
	 let theta = float i in 
	 (f theta) -| (f (theta +. 1.)))

  let nb_binding_regions_by_threshold () = 
    Array.map
      (fun theta -> theta, nbregions (map (binding_in_some_condition theta)))
      [| 1. ; 3. ; 6. ; 9. ; 12. |]

  let nb_confirmed_binding_regions_by_threshold () = 
    Array.map
      (fun theta -> theta, nbregions (map (confirmed_binding_in_some_condition theta)))
      [| 1. ; 3. ; 6. ; 9. ; 12. |]

  let auc motif selected = 
    let selected = filter selected regions in
    if Array.length selected = 0 then nan
    else (
      let bed = Bed.of_locations 
	(Filename.temp_file "guiz" ".bed")
	(Array.enum selected) in
      let peakseq = Ucsc.fasta_of_bed Genome.id bed in 
      let ctrlseq = Fasta.enumerate (Ucsc.control_sequences Genome.id peakseq) in
      Motif_library.(roc_auc (eval_motif motif peakseq ctrlseq)#value)
    )

  let auc_for_simple_binding_regions () =
    Array.map
      (fun theta -> 
	 let regions = map (binding_in_some_condition theta) in
	 theta, nbregions regions, auc balmer_dr0125 regions)
      [| 1. ; 3. ; 6. ; 9. ; 12. |]
      
  let auc_for_confirmed_regions () = 
    Array.map
      (fun theta -> 
	 let regions = map (confirmed_binding_in_some_condition theta) in
	 theta, nbregions regions, auc balmer_dr0125 regions)
      [| 1. ; 3. ; 6. ; 9. ; 12. |]

  let auc_for_simple_binding_regions_by_steps () = 
    Array.map (fun regions -> nbregions regions, auc balmer_dr0125 regions) (etages (fun theta -> map (binding_in_some_condition theta)))

  let auc_for_confirmed_binding_regions_by_steps () = 
    Array.map (fun regions -> nbregions regions, auc balmer_dr0125 regions) (etages (fun theta -> map (confirmed_binding_in_some_condition theta)))

  let auc_for_rar_not_rxr () = 
    Array.map
      (fun theta -> 
	 let regions = map (rar_binding_in_some_condition theta) -| map (rxr_binding_in_some_condition theta) in
	 theta, nbregions regions, auc balmer_dr0125 regions)
      [| 1. ; 3. ; 6. ; 9. ; 12. |]

  let auc_for_rar_not_rxr2 () = 
    Array.map
      (fun theta -> 
	 let regions = map (rar_binding_in_some_condition theta) -| map (rxr_binding_in_some_condition 3.) in
	 theta, nbregions regions, auc balmer_dr0125 regions)
      [| 1. ; 3. ; 6. ; 9. ; 12. |]
    
  let threshold_against_input = 6.
    
  let selected = Array.mapi (fun i _ -> confirmed_binding_in_some_condition threshold_against_input i) regions

  let rar_series = [
    `F9_WT_panRAR_1, `F9_ATRA_panRAR_1 ;
    `F9_ATRA_panRAR_1, `F9_ATRA24_panRAR_1 ;
    `F9_ATRA24_panRAR_1, `F9_ATRA48_panRAR_1 ;
  ]

  let rxr_series = [
    `F9_WT_panRXR_1,     `F9_ATRA2_panRXR_1 ;
    `F9_ATRA2_panRXR_1,  `F9_ATRA24_panRXR_1 ;
    `F9_ATRA24_panRXR_1, `F9_ATRA48_panRXR_1 ;
  ]

  let variation_score x y i = 
    let treatment = (x --> read_counts).(i) 
    and control = (y --> read_counts).(i) 
    and treatment_size = x --> library_size 
    and control_size = y --> library_size in
    let sign = 
      if float treatment /. float treatment_size < float control /. float control_size 
      then 1. else -1.
    and value = pmt_control ~treatment ~control ~treatment_size ~control_size in
    sign *. value

  let raw_variation x y i = 
    let treatment = (x --> read_counts).(i) 
    and control = (y --> read_counts).(i) 
    and treatment_size = x --> library_size 
    and control_size = y --> library_size in
    float treatment /. float treatment_size *. 1e6,
    float control /. float control_size *. 1e6 

  let raw_raw_variation x y i = 
    let treatment = (x --> read_counts).(i) 
    and control = (y --> read_counts).(i) in
    treatment, control

  let raw_raw_variation_profile series i =
    List.map (fun (t1,t2) -> raw_raw_variation t1 t2 i) series


  let discrete theta x = 
    if abs_float x > theta then 
      if x > 0. then `P else `M
    else `Z

  let raw_variation_profile series i =
    List.map (fun (t1,t2) -> raw_variation t1 t2 i) series

  let score_profile series i =
    List.map (fun (t1,t2) -> variation_score t1 t2 i) series

  let variation_profile series theta i =
    List.map (fun (t1,t2) -> discrete theta (variation_score t1 t2 i)) series

  let coherent_profile alpha theta i = 
    variation_profile rar_series alpha i = variation_profile rxr_series theta i

  let coherent_non_null_profile alpha theta i = 
    coherent_profile alpha theta i && variation_profile rar_series alpha i <> [`Z;`Z;`Z]

  let nbcoherent_profiles alpha theta = nbregions (_and selected (map (coherent_profile alpha theta)))

  let nbcoherent_non_null_profiles alpha theta = nbregions (_and selected (map (coherent_non_null_profile alpha theta)))
    
  let grid_search () = 
    Array.init 6 
      (fun i -> Array.init 6
	 (fun j ->
	    nbcoherent_non_null_profiles (float i) (float j)))
		     
  let plot points = 
    let rp = R.make () in
    R.iv rp "x" points ;
    R.c rp "pdf('delme.pdf')" ;
    R.c rp "plot(x)" ;
    R.c rp "dev.off()" ;
    R.close rp ;
    Sle.shell "evince delme.pdf"
  

end


(* CORRELATION ENTRE VALEURS D'EXPRESSION dans deux échantillons *)
let linear_regression x y path = GeneAnnotation.(
  let rp = R.make () in
  let genes = ensembl#value in
  let x = Array.map (fun g -> let v,_,_ = x --> g.cufflinks_expr_level in v) genes
  and y = Array.map (fun g -> let v,_,_ = y --> g.cufflinks_expr_level in v) genes in
  R.v rp "x" x ;
  R.v rp "y" y ;
  R.c rp "pdf('%s')" path ;
  R.c rp "plot(lm(y ~ x))" ;
  R.c rp "dev.off()" ;
  R.close rp
)
