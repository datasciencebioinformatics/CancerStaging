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
# Paired samples                                                                                                         
paired_sample_df<-data.frame(normal=c(),tumor=c(),case=c())                                                              
                                                                                                                         
# For each case, find the pairs                                                                                          
for (case in unique(merged_data_patient_info_data$case))                                                                 
{                                                                                                                        
    # All samples for case id = "case"                                                                                   
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$case==case,]                               
                                                                                                                         
    # Take the tumor samples                                                                                           
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$patient_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$patient_id)
                {
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case),paired_sample_df)
                }
            }                
    }
}
#######################################################################################################################################
# Take control and normal samples
samples_Tumor  <-colData[colData$tumor_normal=="Primary Tumor","patient_id"]
samples_Normal  <-colData[colData$tumor_normal=="Solid Tissue Normal","patient_id"]
#######################################################################################################################################
# folchange=Expr(Stage i)/Expr(Stage ii and II)
# folchange=rowMeans(unstranded_rpkm[,paired_sample_df$tumor])-rowMeans(unstranded_rpkm[,paired_sample_df$normal])
# folchange=rowMeans(unstranded_rpkm[,paired_sample_df$tumor])/rowMeans(unstranded_rpkm[,paired_sample_df$normal])
# log2change=log(folchange,2)	
# Log(FC) = mean(log2(Group1)) / mean(log2(Group2))
# Log(FC) = log2(mean(Group1)/(mean(Group2)))
# Paired t-test, RPKM of paired tumor/normal samples
# Plot with 15208 genes.
# Log2foldchange
LOG_CONSTANT=0.001
log2change_paired=log( (rowMeans(unstranded_rpkm[,paired_sample_df$tumor]+LOG_CONSTANT)/rowMeans(unstranded_rpkm[,paired_sample_df$normal]+LOG_CONSTANT)),2)	
log2change=log( (rowMeans(unstranded_rpkm[,samples_Tumor]+LOG_CONSTANT)/rowMeans(unstranded_rpkm[,samples_Normal]+LOG_CONSTANT)),2)	

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
