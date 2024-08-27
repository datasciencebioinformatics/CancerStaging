####################################################################################################################
# Save normalized data                                                                                             #
unstranded_edgeR_rpkm_file     <-         paste(output_dir,"unstranded_edgeR_rpkm.tsv",sep="")			               #
unstranded_dgelist_rpkm_file   <-         paste(output_dir,"unstranded_dgelist_rpkm.tsv",sep="")		               #
unstranded_NOISeq_rpkm_file    <-         paste(output_dir,"unstranded_NOISeq_rpkm.tsv",sep="")	  		             #
unstranded_NOISeq_TMM_file     <-         paste(output_dir,"unstranded_NOISeq_TMM.tsv",sep="")	 		               #
merged_data_patient_info_file  <-         "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv"     #
###########################################################################################################################
unstranded_edgeR_rpkm_data       <- read.table(file = unstranded_edgeR_rpkm_file, sep = '\t', header = TRUE,fill=TRUE)    #
unstranded_dgelist_rpkm_data     <- read.table(file = unstranded_dgelist_rpkm_file, sep = '\t', header = TRUE,fill=TRUE)  #
unstranded_NOISeq_rpkm_data      <- read.table(file = unstranded_NOISeq_rpkm_file, sep = '\t', header = TRUE,fill=TRUE)   #
unstranded_NOISeq_TMM_data       <- read.table(file = unstranded_NOISeq_TMM_file, sep = '\t', header = TRUE,fill=TRUE)    #
merged_data_patient_info         <- read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE) #
###########################################################################################################################
# Save tumor samples
tumor_samples<-colData[which(colData$tissue_type=="Tumor"),"patient_id"]

# Select only the tumor genes
tumor_genes<-log2change_tumor_control[intersect(which(log2change_tumor_control$FDR<=0.05), which(log2change_tumor_control$log2change>=threshold_tumor)),"gene"]

# Filter by log2folchage
unstranded_rpkm<-unstranded_rpkm[tumor_genes,]

# Filter by RPKM
unstranded_data_filter<-unstranded_rpkm[rowMeans(unstranded_rpkm)>threshold_rpkm,]
###########################################################################################################################
cat(print(paste("\nNumber of up-regulated tumor-genes :",dim(unstranded_data_filter)[1])),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
