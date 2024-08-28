###########################################################################################################################
# Save normalized data                                                                                             
merged_data_patient_info_count<-merged_data_patient_info_count[merged_data_patient_info_count$sample_id %in% colnames(df_reads_count_all_projects_fpkm),]
###########################################################################################################################
# All tumor and control samples
unpaired_tumor_samples  <-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor","sample_id"]
unpaired_control_samples<-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal","sample_id"]

unpaired_tumor_samples   %in%   colnames(df_reads_count_all_projects_fpkm)
unpaired_control_samples %in%   colnames(df_reads_count_all_projects_fpkm)

# Paired samples only
paired_normal_samples    <- paired_sample_df$normal
paired_tumor_samples     <- paired_sample_df$tumor
#######################################################################################################################################
# folchange=Expr(Stage i)/Expr(Stage ii and II)
# Paired t-test, RPKM of paired tumor/normal samples
# Plot with 15208 genes.
# Log2foldchange
LOG_CONSTANT=0.001
log2change       =log( (rowMeans(df_reads_count_all_projects_fpkm[,unpaired_tumor_samples]+LOG_CONSTANT)/rowMeans(df_reads_count_all_projects_fpkm[,unpaired_control_samples]+LOG_CONSTANT)),2)	
log2change_paired=log( (rowMeans(df_reads_count_all_projects_fpkm[,paired_tumor_samples]+LOG_CONSTANT)/rowMeans(df_reads_count_all_projects_fpkm[,paired_normal_samples]+LOG_CONSTANT)),2)	

# log2change data
log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))
log2change_tumor_control_paired=na.omit(data.frame(gene=names(log2change_paired),log2change=log2change_paired))

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Pvalue<-1
log2change_tumor_control_paired$Pvalue<-1

# For each genes in the tabe
for (gene in log2change_tumor_control$gene)
{
	# Take p-value
	log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(df_reads_count_all_projects_fpkm[gene,unpaired_tumor_samples]), y=as.numeric(df_reads_count_all_projects_fpkm[gene,unpaired_control_samples]), paired = FALSE, alternative = "two.sided")$p.value	
	log2change_tumor_control_paired[gene,"Pvalue"]<-t.test(x=as.numeric(df_reads_count_all_projects_fpkm[gene,paired_tumor_samples]), y=as.numeric(df_reads_count_all_projects_fpkm[gene,paired_normal_samples]), paired = TRUE, alternative = "two.sided")$p.value	
}
# FRD 
log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="fdr")
log2change_tumor_control_paired$FDR<-p.adjust(log2change_tumor_control_paired$Pvalue, method="fdr")
#######################################################################################################################################
colnames(log2change_tumor_control)       <- c("gene","log2change_all_samples","pvalue_all_samples","fdr_all_samples")                 #
colnames(log2change_tumor_control_paired)<- c("gene","log2change_paired","pvalue_paired","fdr_paired")                                #
#######################################################################################################################################
logchange_tumor_control<-merge(log2change_tumor_control,log2change_tumor_control_paired,by="gene")                                    #
#######################################################################################################################################








