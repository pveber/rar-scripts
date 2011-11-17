mkdir -p mapability/mm9
cd mapability/mm9
wget ftp://hgdownload.cse.ucsc.edu/gbdb/mm9/bbi/crgMapabilityAlign36mer.bw
bigWigToWig crgMapabilityAlign36mer.bw crgMapabilityAlign36mer.wig
grep -v '^#' crgMapabilityAlign36mer.wig | gawk '{ if( $4 >= 0.99 ) print }' > crgMapabilityAlign36mer.light.wig
cd ../..

mkdir -p conservation/mm9
cd conservation/mm9
wget ftp://hgdownload.cse.ucsc.edu/gbdb/mm9/multiz30way/phastCons30wayPlacental.wib
hgWiggle -db=mm9 phastCons30wayPlacental > phastCons30wayPlacental.wig
cd ../..

