#How to install R for jupyter notebook: https://www.continuum.io/blog/developer/jupyter-and-conda-r
#Particular codes to run under the anaconda command shell: conda install -c r r-essentials
#Then you can "choose" the kernel at the top of jupyter notebook.
#Sleuth homepage: https://pachterlab.github.io/sleuth/about
#How to install Sleuth: https://pachterlab.github.io/sleuth/download.html
#You can always type ?function to learn more about what's going on

#Enter a high memory comp using qrsh -l m_mem_free=32G before entering R.
#Otherwise sleuth won't run and you get an error when using sleuth_prep
#Enter "R" at the terminal

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
library("dplyr")


#'R' syntax:
#<- defines a variable; use "/" for file locations.  You can also just use =
#sample_id is storing a filepathway.  The way this works in R is to use the file.path function where every common denotes a subfolder or ultimately a file.
#kal_dirs creates a variable that holds a list of all the directories that hold kallisto information.
#Metadata denotes where the metadata (i.e. the names of the samples + the condition status to be compared in the Wald analysis) is stored.

base_dir = "C:/Users/philg/Dropbox/notebooks/ipython/human_genomes"
sample_id = file.path(base_dir,"kallisto","SRP070710")
kal_dirs = list.dirs(sample_id, recursive=F) #this stops it from listing SRP070710 as a stand alone directory
metadata = file.path(base_dir,"GSE78220sleuthmeta.txt")

#Create a table by using the read.table or read.csv function to read the metadata.
#Now append a new column to the table that describes the experiment, labeled 'path'. **This ensures that samples can be associated with kallisto quantification

s2c = read.table(file.path(metadata), header = TRUE,sep = "\t", stringsAsFactors=FALSE)
s2c = dplyr::select(s2c,sample,tumor_response)
s2c = dplyr::mutate(s2c,path=kal_dirs)
so = sleuth_prep(s2c, ~condition, num_cores=1, target_mapping = t2g) #specify num_cores is needed for the windows PC version

#Note the addition of the "target_mapping", this maps the genes to the transcript reads.

so = sleuth_fit (so)

#Now find the object that sleuth_prep has created i.e. the "so" variable. Use the sink() function to put this into a file that can then be read.

sink("mymod")
models(so)
sink()
mymodel=scan("mymod",what="list")
condition=mymodel[length(mymodel)]

#from here we find our formula and coefficients
#Take this and perform a "wald test" on *each* coefficients in the format of
# sleuth_wt(obj, "which_beta").  Again store this in "so"

so = sleuth_wt(so,condition)

#You will then extract the wald test results using sleuth_results(obj, test)
#store in table and then order

results_table = sleuth_results(so,condition)
results_ordered=results_table[order(results_table$qval),]

#Create a contingency table of the counts at each combination of factor levels following:
#table(one_or_more_obj_intrepreted_as_factors_or_dataframe)

table(results_ordered$qval<=as.numeric(qi))
#presumably this makesa table only if qval <=stated

write.table(subset(results_ordered,qval<=as.numeric(qi)),file=file.path(sample_id,paste("sleuth.melanoma_transcripts.qval_",qi,"_2condition.txt",sep="")),sep="\t",row.names=F,quote=F)
#This creates a table file to hold these results and tells it where to store it

#now create a gene rankfile as well as a table with the kallisto tpm counts which will be stored in the folder
sleuth_genes = sleuth_gene_table(so, condition, test_type='wt')
sleuth_genes = filter(sleuth_genes, qval<=as.numeric(qi))
sleuth_genes = sleuth_genes[order(sleuth_genes$qval),]
write.table(sleuth_genes,file.path(sample_id,paste("sleuth_genes_",qi,"_2rank.txt", sep="")),sep="\t",row.names=F,quote=F)


tpm_counts=kallisto_table(so)
write.table(tpm_counts,file.path(sample_id,paste("tpm_",qi,"_2counts.txt", sep="")),sep="\t",row.names=F,quote=F)
