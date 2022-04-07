library(facets)

setwd("/tmp")


set.seed(1234)


# datafile = system.file("extdata", "stomach.csv.gz", package="facets")
datafile = system.file("../../../../../../tmp", "data.txt.gz", package="facets")
#head(read.csv(datafile)[,c(1:2,5:12)])




rcmat = readSnpMatrix(datafile)
xx = preProcSample(rcmat,ndepth=35,snp.nbhd=500,gbuild="hg38",cval=25)
#xx = preProcSample(rcmat,ndepth=35,snp.nbhd=250,gbuild="hg38",cval=25)(WES)
oo=procSample(xx,cval=400)
#oo=procSample(xx,cval=150)(WES)
oo$dipLogR

fit=emcncf(oo)

fit$cncf

fit$purity

fit$ploidy

plotSample(x=oo,emfit=fit)

logRlogORspider(oo$out, oo$dipLogR)


