#How to install R for jupyter notebook: https://www.continuum.io/blog/developer/jupyter-and-conda-r
#Particular codes to run under the anaconda command shell: conda install -c r r-essentials
#Then you can "choose" the kernel at the top of jupyter notebook.
#Sleuth homepage: https://pachterlab.github.io/sleuth/about
#How to install Sleuth: https://pachterlab.github.io/sleuth/download.html
#You can always type ?function to learn more about what's going on

#Enter a high memory comp using qrsh -l m_mem_free=32G before entering R.
#Otherwise sleuth won't run and you get an error when using sleuth_prep
#Enter "R" at the terminal
###\\\ HOW TO USE \\\###

#This will evaluate 3 conditions instead of 2 conditions.  The best example would be
#For a gene set where you have CR, PR, Progressive.
#You will need to arrange your metadata appropriately before using this script.

###\\\\\\###

#First time install code:
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#biocLite("devtools")
#devtools::install_github("pachterlab/sleuth")
#install.packages("shiny")
#install.packages("dplyr")

#load sleuth and shiny

library("sleuth")
library("shiny")

userinput = function() {
  val = readline(prompt="Provde the q-value integer you wish to use for the wald test (e.g. 0.05, 0.1) ")
  print(val)
}

qi=userinput()

#'R' syntax:
#<- defines a variable; use "/" for file locations.  You can also just use =
#sample_id is storing a filepathway.  The way this works in R is to use the file.path function where every common denotes a subfolder or ultimately a file.
#kal_dirs creates a variable that holds a list of all?dpl the directories that hold kallisto information.
#Metadata denotes where the metadata (i.e. the names of the samples + the condition status to be compared in the Wald analysis) is stored.

base_dir = "D:/"
sample_id = file.path(base_dir,"kallisto","SRP070710")
kal_dirs = list.dirs(sample_id, recursive=F) #this stops it from listing SRP070710 as a stand alone directory
metadata = file.path(sample_id,"kallisto_meta.csv")

#The code below brings in the "target_id to gene" conversion table from biomaRt and write it into the sample folder.
#If the file was already downloaded, then it will load the file.  Otherwise it downloads a new one.

if (file.exists(file.path(sample_id,"targetid_to_genes.txt"))) {
  t2g = read.table(file.path(sample_id,"targetid_to_genes.txt"), header=TRUE, colClasses="character")
} else {
  mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
  t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
  t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  write.table(t2g, file=file.path(sample_id,"targetid_to_genes.txt"),sep="\t",row.names=F,quote=F)
}


#Create a table by using the read.table or read.csv function to read the metadata.
#Now append a new column to the table that describes the experiment, labeled 'path'. **This ensures that samples can be associated with kallisto quantification

s2c = read.table(file.path(metadata), header = TRUE,sep = ",", stringsAsFactors=FALSE)
s2c = dplyr::select(s2c,sample,tumor_response)
s2c = dplyr::mutate(s2c,path=kal_dirs)
so = sleuth_prep(s2c, ~tumor_response, num_cores=1, target_mapping = t2g) #specify num_cores is needed for the windows PC version

#Note the addition of the "target_mapping", this maps the genes to the transcript reads.
so = sleuth_fit (so)

#Now find the object that sleuth_prep has created i.e. the "so" or "sleuth_object" variable.  Let's automate this:
#use the "sink" function to set stdout to a file called "mymod" then run models(so)
#Run sink() to return stdout default
#use "scan" and define the output as a list (as that is what it is), then set the conditions
#Condition0 should be the bottom of the list, and so-on.  If you are only comparing 2 conditions then you just need condition0

sink("mymod")
models(so)
sink()
mymodel=scan("mymod",what="list")
condition0=mymodel[length(mymodel)]
condition1=mymodel[length(mymodel)-1]
#Now find the object that sleuth_prep has created i.e. the "so" variable

#from here we find our formula and coefficients
#Take this and perform a "wald test" on *each* coefficients in the format of
# sleuth_wt(obj, "which_beta").  Again store this in "so"

so0 = sleuth_wt(so,condition0)
#note that I stored this in a new variable so that I could reuse "so" later
#You will then extract the wald test results using sleuth_results(obj, test)
#store in table and then order

results_table = sleuth_results(so0,condition0)
results_ordered=results_table[order(results_table$qval),]

#Create a contingency table of the counts at each combination of factor levels following:
#table(one_or_more_obj_intrepreted_as_factors_or_dataframe)

table(results_ordered$qval<=as.numeric(qi))
#presumably this makesa table only if qval <=.05

write.table(subset(results_ordered,qval<=as.numeric(qi)),file=file.path(sample_id,paste("sleuth.melanoma_transcripts.qval_",qi,"_3condition0.txt",sep='')),sep="\t",row.names=F,quote=F)
#This creates a table file to hold these results and tells it where to store it

#Now do it again for the next conditions

so1 = sleuth_wt(so,condition1)
results_table = sleuth_results(so1, condition1)
results_ordered=results_table[order(results_table$qval),]
table(results_ordered$qval<=as.numeric(qi))
write.table(subset(results_ordered,qval<=as.numeric(qi)),file=file.path(sample_id,paste("sleuth.melanoma_transcripts.qval_",qi,"_3condition1.txt",sep='')),sep="\t",row.names=F,quote=F)

#now create a gene rankfile as well as a table with the kallisto tpm counts which will be stored in the folder

sleuth_genes_c0 = sleuth_gene_table(so0, condition0, test_type='wt')
sleuth_genes_c1 = sleuth_gene_table(so1, condition1, test_type='wt')
head(sleuth_genes_c0)
library(dplyr)
sleuth_genes_c1 = filter(sleuth_genes_c1, qval<=as.numeric(qi))
sleuth_genes_c0 = filter(sleuth_genes_c0, qval<=as.numeric(qi))
sleuth_genes = rbind(sleuth_genes_c0, sleuth_genes_c1)
sleuth_genes = sleuth_genes[order(sleuth_genes$qval),]
write.table(sleuth_genes,file.path(sample_id,paste("sleuth_genes_",qi,"_3rank.txt", sep="")),sep="\t",row.names=F,quote=F)
tpm_counts=kallisto_table(so)
write.table(tpm_counts,file.path(sample_id,paste("tpm_",qi,"_3counts.txt", sep="")),sep="\t",row.names=F,quote=F)
