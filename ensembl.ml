open Guizmin

module Guizmin_plugin = struct
  let string_of_org = function
  `mouse -> "mouse"
  let transcript_url = function
  `mouse -> "ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.NCBIM37.63.gtf.gz"

  let transcripts x = target 
    [ "mkdir -p $@_download" ;
      sp "cd $@_download && wget %s" (transcript_url x) ;
      "cd $@_download && gunzip *.gz" ;
      "D=${PWD} && cd $@_download && gawk '{print \"chr\"$0}' *.gtf | sed 's/chrMT/chrM/g' | grep -v 'ENSMUST00000127664' > ${D}/$@" ;
      "rm -rf $@_download" ]
    []
end

let promoters gtf = assert false
