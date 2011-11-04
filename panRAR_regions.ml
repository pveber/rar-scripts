let conditions = [ 
  `F9_WT_panRAR_1 ; 
  `F9_ATRA_panRAR_1 ;
  `F9_ATRA24_panRAR_1 ;
  `F9_ATRA48_panRAR_1 ;
  `F9_WT_panRXR_1 ; 
  `F9_ATRA2_panRXR_1 ;
  `F9_ATRA24_panRXR_1 ; 
  `F9_ATRA48_panRXR_1
]

let paired_conditions = [
  `F9_WT_panRAR_1, `F9_WT_panRXR_1 ;  
  `F9_ATRA_panRAR_1, `F9_ATRA2_panRXR_1 ;
  `F9_ATRA24_panRAR_1, `F9_ATRA24_panRXR_1 ;
  `F9_ATRA48_panRAR_1, `F9_ATRA48_panRXR_1
]

let all_of_them = 
  Binding_loci.of_macs_targets
    ~radius:250
    (List.map B.Chipseq.macs_peaks conditions)

