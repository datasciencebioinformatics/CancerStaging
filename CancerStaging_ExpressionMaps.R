# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_scheme<-"tpm"

# Construnction of 3d coordinates.
# X=T2           : Thymine composition in second codon position (T2)
# Y=connectivity : interactome_data                 # Carels checked    .
# Z=expression   : reads_count_all_projects         # Use the fpkm
# Reading the contents of TSV file using read_tsv() method
Interactomes_GC3_T2_file <-"/home/felipe/Documents/github/CancerStaging/Interactomes_GC3_T2.csv"

# Read data
Interactomes_GC3_T2_data <-read.table(file = Interactomes_GC3_T2_file, sep = '\t', header = TRUE,fill=TRUE) 

# Rename second collumn
colnames(Interactomes_GC3_T2_data)[2]<-"SYMBOL"

# Merge Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL
Interactomes_GC3_T2_merged<-merge(geneLength_ENTREZID_ENSEMBL,Interactomes_GC3_T2_data,by="SYMBOL")

# Accross cancer types.
# The hypohtoses is to put the clusters ordered by entropy in the 3d maps                   .
# The hypotheseis is to order the clusters by entropy in the 3d map                         .
# Amont the 2d coordinates this can have low signal, in a three map can have a stroger sinal.
expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

# I need take the intersection of ids from the interactome table and expression table.
# The expression table has unique ids but it seems the Interactomes_GC3_T2_merged is duplicated (check).
# Uniprot entries have multiple ENSEMBL. Criteria is to take UNIRPOTKB with highest conections.
# To Do : take UNIRPOTKB with highest conections.
# Trasnsform data.frame to a data.table
Interactomes_GC3_T2_merged<-as.data.table(Interactomes_GC3_T2_merged)

# Take the occurance with greatest number of connections
Interactomes_GC3_T2_merged<-data.frame(Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged[, .I[which.max(Conections)], by=ENSEMBL]$V1])

# Set ensembl ids
ENSEMBL_ids<-unique(intersect(rownames(expression_table_normalized),Interactomes_GC3_T2_merged$ENSEMBL))

# Interactomes_GC3_T2_merged
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$ENSEMBL %in% ENSEMBL_ids,c("T2","GC3","Conections","ENSEMBL")]

# "Merge expression table" and "Interactomes_GC3_T2_merged"
merged_expression_interactomes<-cbind(expression_table_normalized[ENSEMBL_ids,],Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$ENSEMBL %in% ENSEMBL_ids,c("T2","GC3","Conections","ENSEMBL")])

# Melt data.frame 
melt_expression_interactomes <- melt(merged_expression_interactomes, id=c("T2","GC3","Conections","ENSEMBL"))

# Set name of normalization schme
colnames(melt_expression_interactomes)[6]<-normalization_scheme

# Packages and data use throught 
library(metR)
library(ggplot2)
# FindClusters_resolution
png(filename=paste(output_dir,"geom_contour_filled.png",sep=""), width = 24, height = 24, res=600, units = "cm")
  ggplot(melt_expression_interactomes, aes(T2, GC3, z = tpm)) +
    geom_contour_fill() +
    geom_contour(color = "black", size = 0.1)
dev.off()
