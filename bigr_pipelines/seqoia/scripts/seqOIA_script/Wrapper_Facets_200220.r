

library(facets)

setwd("/tmp")


set.seed(1234)


# datafile = system.file("extdata", "stomach.csv.gz", package="facets")
datafile = system.file("../../../../../../tmp", "data.txt.gz", package="facets")
#head(read.csv(datafile)[,c(1:2,5:12)])



print("rcmat start")
rcmat = readSnpMatrix(datafile)
print("preprocSample")
xx = preProcSample(rcmat,ndepth=35,snp.nbhd=500,gbuild="hg38",cval=25)
#xx = preProcSample(rcmat,ndepth=35,snp.nbhd=250,gbuild="hg38",cval=25)(WES)
print("procSample")
oo=procSample(xx,cval=400)
#oo=procSample(xx,cval=150)(WES)
print("oo")
oo$dipLogR
print("fit")
fit=emcncf(oo)

## Parsing de la table ##


## Altere ou non ##
Altere_ou_non = rep ("Altere",length(fit$cncf$tcn.em))
pos_Normal= which(fit$cncf$tcn.em  == 2  & fit$cncf$lcn.em == 1 )
Altere_ou_non[pos_Normal] = "Normal"
###################

## Quel_alteration ##
Quelle_alteration = rep ("Normal",length(fit$cncf$tcn.em))
#####################

## Gain ##
pos_Gain = which(fit$cncf$tcn.em  > 2 )
Quelle_alteration[pos_Gain] = "Gain"

## Pertes ##
pos_Del = which(fit$cncf$tcn.em  == 1  )
Quelle_alteration[pos_Del] = "Deletion"

pos_DelH = which(fit$cncf$tcn.em  ==  0 )
Quelle_alteration[pos_DelH] = "Del_Homo"

## LOH ##

pos_LoH= which(fit$cncf$tcn.em  == 2  & fit$cncf$lcn.em == 0 )
Quelle_alteration[pos_LoH] = "LoH"
taille_segment = fit$cncf$end - fit$cncf$start 
Table_out = cbind( fit$cncf,length= taille_segment , Altere_ou_non, Quelle_alteration, Purity  = fit$purity , Infered_ploidy = fit$ploidy )

reorder_colnames = c("chrom","start","end","length","num.mark","nhet","cnlr.median","mafR","segclust",
                     "cnlr.median.clust","mafR.clust",
                    "cf.em","tcn.em","lcn.em","Altere_ou_non","Purity","Infered_ploidy")

Table_out_reordered = Table_out[,reorder_colnames]



write.table(Table_out_reordered,"Output_Facets.bed",sep = "\t",row.names=FALSE)
  
plotSample(x=oo,emfit=fit)
logRlogORspider(oo$out, oo$dipLogR)
