source("https://bioconductor.org/biocLite.R")
biocLite("sva")
browseVignettes("sva")
#biocLite("bladderbatch")
#biocLite("pamr")
library(sva)
#library(bladderbatch)
#data(bladderdata)
library(pamr)
library(limma)
library(dplyr)

##Example from vignette##
#pheno = pData(bladderEset)

#edata = exprs(bladderEset)

#mod = model.matrix(~as.factor(cancer), data=pheno)
#mod0 = model.matrix(~1, data=pheno)
#n.sv = num.sv(edata, mod, method="leek")
#n.sv
#svobj = sva(edata, mod, mod0, n.sv=n.sv)
#batch = pheno$batch
#modcombat = model.matrix(~1, data=pheno)
#combat_edata = ComBat(dat=edata, batch=batch, mod=mod0,par.prior=TRUE, prior.plots = FALSE)

#pValuesComBat = f.pvalue(combat_edata,mod,mod0)
#qValuesComBat = p.adjust(pValuesComBat, method="BH")

###Now My Code###
setwd("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes")
pheno = read.csv("gsemergeattrib.csv", row.names = 1, header = TRUE, check.names = FALSE) #checknames false prevents it from replacing "-" with "."
rownames(pheno) = pheno$sample_id
pheno = dplyr::select(pheno,response,batch) #selects only the latter two columns from pheno
edata = read.csv("gsemerge_edata.csv", row.names = 1, header = TRUE, check.names = FALSE)
rownames(edata) = edata$gene_id
edata = dplyr::select(edata, -gene_id)
mod = model.matrix(~as.factor(response), data=pheno)
mod0 = model.matrix(~1, data=pheno)
n.sv = num.sv(edata, mod, method = "leek")
n.sv
svobj = sva(as.matrix(edata), mod, n.sv=n.sv)
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=as.matrix(edata), batch=batch, mod=mod0, par.prior=TRUE, prior.plots = FALSE)

pvalues = f.pvalue(combat_edata,mod,mod0)
qvalues = p.adjust(pvalues, method="BH")
write.csv(as.table(qvalues),"gsemerge_qvalues.csv")
write.csv(as.table(pvalues),"gsemerge_pvalues.csv")
