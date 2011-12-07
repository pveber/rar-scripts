library(MASS)

timepoints <- c(0,2,4,6)
rar_columns <- c('coverage.F9.WT.panRAR.1','coverage.F9.ATRA.panRAR.1',
	         'coverage.F9.ATRA24.panRAR.1','coverage.F9.ATRA48.panRAR.1')

data <- read.table('results/chipseq/annotations/PanRAR_regions_chipseq_annotation.tsv',header=TRUE,sep='\t')


points <- t(apply(data[,rar_columns],1,function(x) { x / sqrt(sum(x * x))}))

wss <- (nrow(points)-1)*sum(apply(points,2,var))
for (i in 2:35) wss[i] <- sum(kmeans(points,centers=i)$withinss)
pdf('results/chipseq/clustering_maison/choix_nombre_de_clusters.pdf')
plot(1:35, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") 
dev.off()


km <- kmeans(points,30)
clustering <- km$cluster

x <- seq(0,6,by=0.1)
interp_points <- apply(points,1,function(y) { return(approx(timepoints,y,xout=x)$y) })

profile <- function(sel) { 
  xpoints <- rep(x,dim(points[sel,])[1])
  ypoints <- as.vector(interp_points[,sel])
  f <- kde2d(xpoints,ypoints)
  filled.contour(f,ylim=c(0, 2 * quantile(ypoints,0.8)),
                 color.palette=colorRampPalette(c("white", "cornflowerblue"),space = "rgb"),
		 nlevels=30,
		 plot.axes={
#		   lines(timepoints,c(median(points[sel,1]),median(points[sel,2]),median(points[sel,3]),median(points[sel,4])),col='#FF0000')
		   boxplot(points[sel,1],at=timepoints[1],add=TRUE)
		   boxplot(points[sel,2],at=timepoints[2],add=TRUE)
		   boxplot(points[sel,3],at=timepoints[3],add=TRUE)
		   boxplot(points[sel,4],at=timepoints[4],add=TRUE)
		 })
  
}

pdf('results/chipseq/clustering_maison/cluster_profiles_30.pdf')
for(i in 1:30) {
  profile(clustering == i)
}
dev.off()


output_tab <- data.frame(
	   chrom = data$chrom,
	   start = data$start,
	   end   = data$end,
	   cl5   = kmeans(points,5)$cluster,
	   cl10  = kmeans(points,10)$cluster,
	   cl15  = kmeans(points,15)$cluster,
	   cl20  = kmeans(points,20)$cluster,
	   cl25  = kmeans(points,25)$cluster,
	   cl30  = kmeans(points,30)$cluster,
	   cl35  = kmeans(points,35)$cluster
)

write.table(output_tab,"results/chipseq/clustering_maison/PanRAR_regions_chipseq_clustering.tsv",quote=FALSE,sep='\t')
