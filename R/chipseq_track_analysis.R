data <- read.table('results/chipseq/annotations/PanRAR_regions_chipseq_annotation.tsv',header=TRUE,sep='\t')
tracks <- read.table('results/chipseq/annotations/PanRAR_regions_chipseq_es.tsv',header=TRUE,sep='\t')

prop <- function(col,sel) {
  sum(col[sel]) / sum(sel)
}

profile <- function(sel,title) { 
    p <- c(
      prop(tracks$Wei.Nanog,sel),
      prop(tracks$Wei.Oct4,sel),
      prop(tracks$Wei.Sox2,sel),
      prop(tracks$Wei.Smad1,sel),
      prop(tracks$Wei.E2f1,sel),
      prop(tracks$Wei.Tcfcp2l1,sel),
      prop(tracks$Wei.Zfx,sel),
      prop(tracks$Wei.Stat3,sel),
      prop(tracks$Wei.Klf4,sel),
      prop(tracks$Wei.c.Myc,sel),
      prop(tracks$Wei.n.Myc,sel),
      prop(tracks$Wei.p300,sel),
      prop(tracks$Wei.Suz12,sel),
      prop(tracks$Wei.Esrrb,sel),
      prop(tracks$Wei.CTCF,sel),
      prop(tracks$Wei.Nanog & tracks$Wei.Oct4 & tracks$Wei.Sox2 & tracks$Wei.Smad1 & tracks$Wei.Stat3,sel),
      prop(tracks$Wei.n.Myc & tracks$Wei.c.Myc & tracks$Wei.E2f1 & tracks$Wei.Zfx,sel)
)

    names <- c(colnames(tracks)[1:15],"MTL1","MTL2")

    barplot(p,main=paste(title,sum(sel),"regions"),names.arg=names,cex.names=0.5,las=2)
}


pdf("chipseq_es_tracks.pdf")
profile(1:13860 > 0,"All")
profile(data$pvalue.F9.WT.panRAR.1 > 9,"Bound in WT")
profile(data$pvalue.F9.WT.panRAR.1 < 0.5, "Not bound in WT")
profile(data$pvalue.F9.WT.panRAR.1 < 0.5 & data$pvalue.F9.ATRA.panRAR.1 > 9, "Strongly induced at 2h")
profile(data$pvalue.F9.WT.panRAR.1 < 0.5 & data$pvalue.F9.ATRA.panRAR.1 < 0.5 & data$pvalue.F9.ATRA24.panRAR.1 > 9, "Bound only after 24h")
profile(data$pvalue.F9.WT.panRAR.1 > 9 & data$pvalue.F9.ATRA.panRAR.1 < 2., "Strongly repressed at 2h")
profile(data$pvalue.F9.WT.panRAR.1 > 9 & data$pvalue.F9.ATRA.panRAR.1 > 9 & data$pvalue.F9.ATRA24.panRAR.1 < 0.5, "Strongly repressed at 24h")

dev.off()
