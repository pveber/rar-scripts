library(gplots)

group_color = function(n) { 
  if(n %% 2 == 0) { 
    return('#000000') 
  } 
  else { 
    return('#BBBBBB') 
  } 
}

plot <- function(repr,tf,meth,k) {
  outputfn <- paste('results/chipseq/justin_clustering',
                    repr,
		    paste('chipseq-',tf,'-',meth,'.',k,'.pdf',sep=''),
		    sep='/')
  datafn <- paste('resources/chipseq/justin_clustering',
                  repr,
                  paste('chipseq-',tf,'-',meth,sep=''),
		  sep='/')
		  
  pdf(outputfn)
  i <- k + 6
  data <- read.table(datafn,header = TRUE)
  rowside <- unlist(lapply(data[[i]][order(data[[i]])],group_color))
  heatmap(as.matrix(data[order(data[[i]]),4:7]),
          NA,NA,
          main=paste(repr,tf,meth,'k = ', k),
          col=redgreen(75),
          RowSideColors=rowside)
  dev.off()
}

dir <- function(dirpath,outdirpath) {
  files <- c('chipseq-rar-cos')
  for(i in 1:length(files)) {

    plot(paste(dirpath,files[i],sep='/'), 
         paste(outdirpath,'/',files[i],'.pdf',sep=''),
	 files[i])
  }
}

all <- function() {
  tfs <- c('rar','rxr')
  meths <- c('cos','hsic')
  reprs <- c('clustering','cl-p-norm','cl-p-value','cl-counts-norm','cl-counts')
  ks <- list(clustering = 2:30,
	     cl_counts = c(10,15,19,24),
             cl_counts_norm = c(10,20,30),
	     cl_p_value = c(10,20,30),
	     cl_p_norm = c(10,20,30))

  for(r in 1:length(reprs)) {
    system(paste('mkdir -p results/chipseq/justin_clustering',reprs[r],sep='/'))
    for(tf in tfs) 
      for(meth in c('cos','hsic','spectral'))
        for(k in ks[[r]]) {
	  plot(reprs[r],tf,meth,k)
        }
  }
}

all()
