library(MASS)

timepoints <- c(0,2,24,48)
rar_columns <- c('coverage.F9.WT.panRAR.1','coverage.F9.ATRA.panRAR.1','coverage.F9.ATRA24.panRAR.1','coverage.F9.ATRA48.panRAR.1')

data <- read.table('results/chipseq/annotations/original_PanRAR_regions_chipseq_annotation.tsv',header=TRUE,sep='\t')
points <- data[,rar_columns]

x <- seq(0,48,by=2)
interp_points <- apply(points,1,function(y) { return(approx(timepoints,y,xout=x)$y) })

profile <- function(repr,tf,meth,k) {
  dir <- paste('results/chipseq/justin_clustering/profiles',repr,sep='/')
  outputfn <- paste(dir,
  		    paste('chipseq-',tf,'-',meth,'.',k,'.pdf',sep=''),
		    sep='/')
  datafn <- paste('resources/chipseq/justin_clustering',
                  repr,
                  paste('chipseq-',tf,'-',meth,sep=''),
		  sep='/')
	 
  system(paste('mkdir -p ',dir))
  clustering <- read.table(datafn,header=TRUE,sep='\t')
  pdf(outputfn)
  for(i in 1:k) {
    sel <- clustering$cl30 == (i - 1)
    xpoints <- rep(x,dim(points[sel,])[1])
    ypoints <- as.vector(interp_points[,sel])
    f <- kde2d(xpoints,ypoints)
    filled.contour(f,ylim=c(0, 2 * quantile(ypoints,0.8)),
	           col=gray.colors(50),nlevels=40,
		   plot.axes={
			   lines(timepoints,c(median(points[sel,1]),
                                              median(points[sel,2]),
					      median(points[sel,3]),
					      median(points[sel,4])),
			         col='#FF0000')
		           boxplot(points[sel,1],at=timepoints[1],add=TRUE)
			   boxplot(points[sel,2],at=timepoints[2],add=TRUE)
		           boxplot(points[sel,3],at=timepoints[3],add=TRUE)
			   boxplot(points[sel,4],at=timepoints[4],add=TRUE)
                   },
	           main=paste('cluster',i,'contains',sum(sel),'regions'))
  }
  dev.off()
}

profile('clustering','rar','cos',30)
profile('clustering','rxr','cos',30)
profile('cl-counts','rar','cos',30)
profile('cl-counts','rxr','cos',30)
profile('cl-p-norm','rar','cos',30)
profile('cl-p-value','rar','cos',30)
profile('cl-counts-norm','rar','cos',30)
