clustering <- read.table("resources/chipseq/justin_clustering/cl-counts/chipseq-rar-cos",header = TRUE)
conservation <- read.table("results/chipseq/annotations/original_PanRAR_regions_conservation_annotation.tsv")

pdf('results/chipseq/conservation_in_jclusters.pdf')
boxplot(conservation$V4 ~ clustering$cl30, main='Conservation in clusters',xlab='cluster membership',ylab='Cumulated phastCons score in region')
dev.off()
