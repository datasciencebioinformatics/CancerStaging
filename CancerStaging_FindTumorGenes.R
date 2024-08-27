###########################################################################################################################
# Save normalized data                                                                                             
unstranded_edgeR_rpkm       <- unstranded_rpkm
unstranded_dgelist_rpkm     <- unstranded_dgelist_rpkm
unstranded_NOISeq_rpkm_data <- unstranded_NOISeq_rpkm
unstranded_NOISeq_TMM       <- unstranded_NOISeq_TMM

unstranded_edgeR_rpkm         <-unstranded_edgeR_rpkm[colnames(unstranded_edgeR_rpkm) %in% merged_data_patient_info_count$sample_id,]
merged_data_patient_info_count<-merged_data_patient_info_count[merged_data_patient_info_count$sample_id %in% colnames(unstranded_edgeR_rpkm),]
###########################################################################################################################
# All tumor and control samples
unpaired_tumor_samples  <-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor","sample_id"]
unpaired_control_samples<-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal","sample_id"]

unpaired_tumor_samples   %in%   colnames(unstranded_edgeR_rpkm)
unpaired_control_samples %in%   colnames(unstranded_edgeR_rpkm)

# Paired samples only
paired_normal_samples    <- paired_sample_df$normal
paired_tumor_samples     <- paired_sample_df$tumor
#######################################################################################################################################
# folchange=Expr(Stage i)/Expr(Stage ii and II)
# Paired t-test, RPKM of paired tumor/normal samples
# Plot with 15208 genes.
# Log2foldchange
LOG_CONSTANT=0.001
log2change       =log( (rowMeans(unstranded_edgeR_rpkm[,unpaired_tumor_samples]+LOG_CONSTANT)/rowMeans(unstranded_edgeR_rpkm[,unpaired_control_samples]+LOG_CONSTANT)),2)	
log2change_paired=log( (rowMeans(unstranded_edgeR_rpkm_data[,samples_Tumor]+LOG_CONSTANT)/rowMeans(unstranded_edgeR_rpkm_data[,samples_Normal]+LOG_CONSTANT)),2)	

# log2change data
log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))
log2change_tumor_control_paired=na.omit(data.frame(gene=names(log2change_paired),log2change=log2change_paired))

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Category<-"insignificant"
log2change_tumor_control_paired$Category<-"insignificant"

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Pvalue<-1
log2change_tumor_control_paired$Pvalue<-1

# For each genes in the tabe
for (gene in log2change_tumor_control$gene)
{
	# Take p-value
	log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_rpkm[gene,samples_Tumor]), y=as.numeric(unstranded_rpkm[gene,samples_Normal]), paired = FALSE, alternative = "two.sided")$p.value	
	log2change_tumor_control_paired[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_rpkm[gene,paired_sample_df$tumor]), y=as.numeric(unstranded_rpkm[gene,paired_sample_df$normal]), paired = FALSE, alternative = "two.sided")$p.value	
}
# FRD 
log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="fdr")
log2change_tumor_control_paired$FDR<-p.adjust(log2change_tumor_control_paired$Pvalue, method="fdr")

# Categorize genes if log2foldchange >= threshold_tumor
log2change_tumor_control[intersect(which(log2change_tumor_control$FDR<=threshold_FDR), which(log2change_tumor_control$log2change>=threshold_tumor)),"Category"]<-paste("Tumor genes", sep="")
log2change_tumor_control_paired[intersect(which(log2change_tumor_control_paired$FDR<=threshold_FDR), which(log2change_tumor_control_paired$log2change>=threshold_tumor)),"Category"]<-paste("Tumor genes", sep="")
#######################################################################################################################################
