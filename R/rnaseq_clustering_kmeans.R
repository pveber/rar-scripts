library(MASS)
library(gplots)

data <- read.table('results/rnaseq/ensembl_rnaseq_annotation.tsv',header=TRUE,sep='\t')

timepoints <- c(0,2,4,6,8,10)
relevant_columns <- c('t0.ATRA.t6.baseMeanA',  't0.ATRA.t6.baseMeanB', 
		      't0.ATRA.t12.baseMeanB', 't0.ATRA.t24.baseMeanB',
		      't0.ATRA.t36.baseMeanB', 't0.ATRA.t48.baseMeanB')

under_threshold <- function(th,x) !(is.na(x) | x > th)

#Differentially expressed filter
de_filter <- under_threshold(1e-3,data$t0.ATRA.t6.padj) | under_threshold(1e-3,data$t0.ATRA.t12.padj) | under_threshold(1e-3,data$t0.ATRA.t24.padj) | under_threshold(1e-3,data$t0.ATRA.t36.padj) | under_threshold(1e-3,data$t0.ATRA.t48.padj)

points <- t(apply(data[de_filter,relevant_columns],1,function(x) { x / sqrt(sum(x * x))}))

wss <- (nrow(points)-1)*sum(apply(points,2,var))
for (i in 2:35) wss[i] <- sum(kmeans(points,centers=i)$withinss)
pdf('results/rnaseq/clustering_maison/intra_cluster_variance.pdf')
plot(1:35, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") 
dev.off()

nbclusters = 25
km <- kmeans(points,nbclusters)
clustering <- km$cluster

x <- seq(0,10,by=0.1)
interp_points <- apply(points,1,function(y) { return(approx(timepoints,y,xout=x)$y) })

profile <- function(sel) { 
  xpoints <- rep(x,dim(points[sel,])[1])
  ypoints <- as.vector(interp_points[,sel])
  f <- kde2d(xpoints,ypoints)
  filled.contour(f,ylim=c(0, 2 * quantile(ypoints,0.8)),
                 color.palette=colorRampPalette(c(143, "red"),space = "rgb"),
		 nlevels=20,
		 plot.axes={
#		   lines(timepoints,c(median(points[sel,1]),median(points[sel,2]),median(points[sel,3]),median(points[sel,4])),col='#FF0000')
		   boxplot(points[sel,1],at=timepoints[1],add=TRUE)
		   boxplot(points[sel,2],at=timepoints[2],add=TRUE)
		   boxplot(points[sel,3],at=timepoints[3],add=TRUE)
		   boxplot(points[sel,4],at=timepoints[4],add=TRUE)
		   boxplot(points[sel,5],at=timepoints[5],add=TRUE)
		   boxplot(points[sel,6],at=timepoints[6],add=TRUE)
		 })
  
}
profile(clustering == 3)

pdf('results/rnaseq/clustering_maison/cluster_profiles_30.pdf')
for(i in 1:nbclusters) {
  profile(clustering == i)
}
dev.off()


output_tab <- data.frame(
	   "gene id" = data$gene.id[de_filter],
	   cl5   = kmeans(points,5)$cluster,
	   cl10  = kmeans(points,10)$cluster,
	   cl15  = kmeans(points,15)$cluster,
	   cl20  = kmeans(points,20)$cluster,
	   cl25  = kmeans(points,25)$cluster,
	   cl30  = kmeans(points,30)$cluster,
	   cl35  = kmeans(points,35)$cluster
)

write.table(output_tab,"results/rnaseq/clustering_maison/PanRAR_regions_chipseq_clustering.tsv",quote=FALSE,sep='\t')

