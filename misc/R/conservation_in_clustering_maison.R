clustering <- read.table("results/chipseq/clustering_maison/PanRAR_regions_chipseq_clustering.tsv",header = TRUE,sep='\t')
conservation <- read.table("results/chipseq/annotations/PanRAR_regions_conservation_annotation.tsv")

pdf('results/chipseq/conservation_in_clustering_maison.pdf')
boxplot(conservation$V4 ~ clustering$cl30, main='Conservation in clusters',xlab='cluster membership',ylab='Cumulated phastCons score in region')
dev.off()
