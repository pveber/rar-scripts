data <- read.table('resources/chipseq/regions/PanRAR_allregions_pmt.tsv',header=TRUE,sep='\t')

rar_max_pvalue <- pmax(data$pvalue.F9.WT.panRAR.1,data$pvalue.F9.ATRA.panRAR.1,data$pvalue.F9.ATRA24.panRAR.1,data$pvalue.F9.ATRA48.panRAR.1)
rxr_max_pvalue <- pmax(data$pvalue.F9.WT.panRXR.1,data$pvalue.F9.ATRA2.panRXR.1,data$pvalue.F9.ATRA24.panRXR.1,data$pvalue.F9.ATRA48.panRXR.1)

y = c()
for(i in 1:15) {
      y[i] = sum((rar_max_pvalue > i) & (rxr_max_pvalue > i))
}
pdf('threshold_choice_inter.pdf')
plot(1:15,y)
dev.off()

y = c()
for(i in 1:15) {
      y[i] = sum((rar_max_pvalue > i) & (rxr_max_pvalue > i)) / sum(rar_max_pvalue > i)
}
pdf('threshold_choice_prop.pdf')
plot(1:15,y)
dev.off()
