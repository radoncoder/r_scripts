#Takes p-values from python mann-whitney and turns them into Q values with file output.

source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library("qvalue")
install.packages("dplyr")
library("dplyr")
gse1 = read.csv(file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes",
                   "GSE93157pvalue.csv"), header = TRUE, sep = ",")
pvalues = gse1["p.value"]
qobj = qvalue(p=pvalues)
qvalues = qobj["qvalues"]
q=qvalues$qvalues$p.value
gse1=dplyr::mutate(gse1,qvalues=q)
write.csv(gse1,file = file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes","GSE93157p-qvalue.csv"), row.names = FALSE)

gse2 = read.csv(file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes",
                          "GSE67501pvalue.csv"), header = TRUE, sep = ",")
pvalues = gse2["p.value"]
qobj = qvalue(p=pvalues)
qvalues = qobj["qvalues"]
q=qvalues$qvalues$p.value
gse2=dplyr::mutate(gse2,qvalues=q)
write.csv(gse2,file = file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes","GSE67501p-qvalue.csv"), row.names = FALSE)

gse3 = read.csv(file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes",
                          "GSE78220pvalue.csv"), header = TRUE, sep = ",")
pvalues = gse3["p.value"]
qobj = qvalue(p=pvalues)
qvalues = qobj["qvalues"]
q=qvalues$qvalues$p.value
gse3=dplyr::mutate(gse3,qvalues=q)
write.csv(gse3,file = file.path("C:/Users/philg/Dropbox/notebooks/ipython/human_genomes","GSE78220p-qvalue.csv"), row.names = FALSE)
