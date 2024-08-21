#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
# To construct sub-interactome networks, pairwise combinations of stage-specific genes are created and then filtered to keep edges overlapping the interctome.
#######################################################################################################################################
# Path to input files
# Interactome file
interactome_file<-"/home/felipe/Documents/github/CancerStaging/Full_Interactome_Flavia.txt"

# EnsemblToUniprotKBconversionList file
EnsemblToUniprotKBconversionList_file<-"/home/felipe/Documents/github/CancerStaging/EnsemblToUniprotKBconversionList.txt"
#######################################################################################################################################
# Read input table
# Gene table
interactome_data <-read.table(file = interactome_file, sep = '\t', header = FALSE,fill=TRUE)         

# Gene EnsemblToUniprotKBconversionList_data
EnsemblToUniprotKBconversionList_data <-read.table(file = EnsemblToUniprotKBconversionList_file, sep = '\t', header = TRUE,fill=TRUE)         

# Rename collumns 
colnames(interactome_data)<-c("Gene1","Gene2")

#dim(interactome_data)
#length(unique(c(interactome_data$Gene1,interactome_data$Gene2)))
#######################################################################################################################################
# Before converstion 
interactome_data_raw<-interactome_data

# Filter tables to keep only the gene entries that are listed in the EnsemblToUniprotKBconversionList
interactome_data<-interactome_data[interactome_data$Gene1 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]
interactome_data<-interactome_data[interactome_data$Gene2 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]

# Create a table for id conversion gene_id and gene_symbol for the genes in the interactome data
gene1_conversion<-merge(interactome_data,EnsemblToUniprotKBconversionList_data,by.x="Gene1", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)
gene_conversion<-merge(gene1_conversion,EnsemblToUniprotKBconversionList_data,by.x="Gene2", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)

# Keep only the collumns of interest- 
# interactome_data : interactome with converted ids
interactome_data<-unique(gene_conversion[,3:4])

# Rename interactome_data collumns
colnames(interactome_data)<-c("Gene1","Gene2")
colnames(interactome_data)<-c("Gene1","Gene2")

#dim(interactome_data)[1]/dim(interactome_data_raw)[1]
#length(unique(c(interactome_data$Gene1,interactome_data$Gene2)))/length(unique(c(interactome_data_raw$Gene1,interactome_data_raw$Gene2)))
dim(interactome_data)[1]
dim(interactome_data_raw)[1]
length(unique(c(interactome_data$Gene1,interactome_data$Gene2)))
length(unique(c(interactome_data_raw$Gene1,interactome_data_raw$Gene2)))

# The IntAct interactome was obtained from the intact-micluster.txt file (version updated December 2017) accessed on January 11, 2018, with 152280 interactions among 15651 gene symbols. After converting gene symbols to ENSEMBL identifiers with EnsemblToUniprotKBconversionList.txt, 148169 interactions (97.3%) and 14492 genes (92.6%) were kept. To calculate the connectivity per gene, we counted the number of times each gene appeared in the interactome. 


# To calculate the connectivity of each, we counted the number of times each gene appeared in the interactome. The gene symbols were translated to ensemble identifiers with EnsemblToUniprotKBconversionList.txt.
connectivity<-table(c(gene_conversion$ENSG.x,gene_conversion$ENSG.y))
