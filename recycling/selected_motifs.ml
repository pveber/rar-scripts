open Batteries
open Genome
open Oregon
open Motif_library
open Printf

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

let balmer_dtrn i j =
  PSSM (`sequence [ balmer_hexamer ; `gap (i,i) ; balmer_hexamer ; `gap (j,j) ; balmer_hexamer ],
	sprintf "BalmerDTR(%d,%d)" i j)

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
    

let sox2_oct4_composite_by_eye = 
  PSSM (`sequence [
	  `base (-1.,0.2,0.2,0.) ;
	  `base (0.3,-1.,-1.,0.3) ;
	  `base (-1.,-1.,-1.,2.) ;
	  `base (-1.,-1.,-1.,2.) ;
	  `base (0.05,0.1,0.7,-1.) ;
	  `base (0.2,-1.,-1.,1.) ;
	  `base (-0.1,0.1,-0.1,0.1) ;
	  `base (1.3,0.1,-1.,0.1) ;
	  `base (-1.,-1.,-1.,2.) ;
	  `base (-1.,-1.,1.3,0.1) ;
	  `base (-1.,0.3,0.1,0.1) ;
	  `base (0.6,-1.,-1.,0.1) ;
	  `base (0.4,-1.,0.1,0.05) ;
	  `base (0.6,0.1,-1.,0.1) ;
	  `base (0.1,0.1,-1.,0.4) ;
	],
	"ugly_sox2_oct4")

let value = Target.V.make
  (object
     method id = "Selected_motif.value[r1]"
     method deps = []
     method build = Array.concat [
       [| sox2_oct4_composite_by_eye |];
       Array.init 6 balmer_drn ;
       Array.init 6 balmer_ern ;
       Array.init 6 balmer_irn ;
       [| balmer_dr125 ; balmer_dr0125 |]
     ]
   end)

let drerir = Array.concat [
  Array.init 10 balmer_drn ;
  Array.init 10 balmer_ern ;
  Array.init 10 balmer_irn ;
]

let trimer = 
  let a = Array.init 10 identity  in
  Array.map
    (fun (i,j) -> balmer_dtrn i j)
    (Sle.SleArray.product a a)

let jaspar_tandems () = Motif_library.(
  motif_tandems ~mingap:0 ~maxgap:200 balmer_dr0125 (of_jaspar_collection Jaspar.core)
)
