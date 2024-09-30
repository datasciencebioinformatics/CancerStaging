#######################################################################################################################
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")                               #
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")                                   #
#######################################################################################################################
# saveRDS                                                                                                             #
saveRDS(object = df_reads_count_all_projects, file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))      #
saveRDS(object = geneLength_ENTREZID_ENSEMBL, file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))      #
saveRDS(object = Interactomes_GC3_T2_merged, file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))        #
                                                                                                                      #
# Restore the object                                                                                                  #
geneLength_ENTREZID_ENSEMBL<-readRDS(file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))               #
Interactomes_GC3_T2_merged<-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))                 #
df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))               #
#######################################################################################################################
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

# Read the Interactomes_GC3_T2_data
# First  collumn : UniprotKB
# Second collumn : SYMBOL
# Third  collumn : GC3
# Fourth collumn : T2
# Fifth  collumn : Conections
Interactomes_GC3_T2_data <-read.table(file = Interactomes_GC3_T2_file, sep = '\t', header = TRUE,fill=TRUE) 

# Set the second collum from GeneSymbol to SYMBOL
colnames(Interactomes_GC3_T2_data)[2]<-"SYMBOL"

# Merge Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL by gene SYMBOL.
# geneLength_ENTREZID_ENSEMBL contains the following collumn:
# ENTREZID   entrezid identifier
# geneLength genelenght calculated by Carels
# SYMBOL     gene symbol
# ENSEMBL    ENSEMBL symbol
Interactomes_GC3_T2_merged<-merge(geneLength_ENTREZID_ENSEMBL,Interactomes_GC3_T2_data,by="SYMBOL")

# Take also expression data from the normalization scheme set by "normalization_scheme"
expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

# Transform the Interactomes_GC3_T2_merged to as.data.table
Interactomes_GC3_T2_merged<-as.data.table(Interactomes_GC3_T2_merged)

# In the file "Interactomes_GC3_T2" each genes can have multiple entries with multiple values for the variable connections.
# Only the occurance with greatest number of Conections will be used.
# Interactomes_GC3_T2_merged
Interactomes_GC3_T2_merged<-data.frame(Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged[, .I[which.max(Conections)], by=ENSEMBL]$V1])

# Set ensembl ids
ENSEMBL_ids<-unique(intersect(rownames(expression_table_normalized),Interactomes_GC3_T2_merged$ENSEMBL))

# Selecte collumns to be extracted, Interactomes_GC3_T2_merged
# "T2","GC3","Conections" and "ENSEMBL"
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$ENSEMBL %in% ENSEMBL_ids,c("T2","GC3","Conections","ENSEMBL")]
####################################################################################################################################################
# Set AveExp to zero 
Interactomes_GC3_T2_merged$AveExp<-0

# Calculate the average expression for the epression of each g
Interactomes_GC3_T2_merged[ENSEMBL_ids,"AveExp"]<-rowMeans(expression_table_normalized[ENSEMBL_ids,])

# merged_expression_interactomes
merged_expression_interactomes<-cbind(expression_table_normalized[ENSEMBL_ids,],Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$ENSEMBL %in% ENSEMBL_ids,c("T2","GC3","Conections","ENSEMBL")])

# Melt data.frame 
melt_expression_interactomes <- melt(merged_expression_interactomes, id=c("T2","GC3","Conections","ENSEMBL"))

# Set colnames
colnames(melt_expression_interactomes)[6]<-normalization_scheme
####################################################################################################################################################
# I have two tables to be use.                                                                                                                     #
# First , a table with with "T2", "GC3", "Conections", "ENSEMBL" "AveExp" per gene                                                                 #
# Interactomes_GC3_T2_merged                                                                                                                       #
# Second, a table with with "T2", "GC3", "Conections", "ENSEMBL" "Exp"    per patient                                                              #
# melt_expression_interactomes                                                                                                                     #
####################################################################################################################################################
Interactomes_GC3_T2_selected                       <-melt_expression_interactomes[,c("T2","GC3",normalization_scheme,"Conections")]
Interactomes_GC3_T2_selected                       <-melt_expression_interactomes[!is.na(melt_expression_interactomes$Conections),]
Interactomes_GC3_T2_selected                       <-melt_expression_interactomes[!is.na(melt_expression_interactomes$T2),]
Interactomes_GC3_T2_selected                       <-melt_expression_interactomes[!is.na(melt_expression_interactomes$GC3),]
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_selected[Interactomes_GC3_T2_selected[,normalization_scheme]>0,]

# FindClusters_resolution
png(filename=paste(output_dir,"geom_contour_melt.png",sep=""), width = 24, height = 24, res=600, units = "cm")  
        scatterplot3d(Interactomes_GC3_T2_selected[,c("T2","GC3",normalization_scheme)], pch = 16)
dev.off()
####################################################################################################################################################
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_merged[,c("T2","GC3","AveExp","Conections")]
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_merged[!is.na(melt_expression_interactomes$Conections),]
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_merged[!is.na(melt_expression_interactomes$T2),]
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_merged[!is.na(melt_expression_interactomes$GC3),]
Interactomes_GC3_T2_selected                       <-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$AveExp>0,]

# FindClusters_resolution
png(filename=paste(output_dir,"geom_contour_merged.png",sep=""), width = 24, height = 24, res=600, units = "cm")  
        scatterplot3d(Interactomes_GC3_T2_merged[,c("T2","GC3","AveExp")], pch = 16)
dev.off()
####################################################################################################################################################
# Plost histogram of T2, GC3 and AveExp 
# Plost histogram of T2, GC3 and TPM
####################################################################################################################################################
