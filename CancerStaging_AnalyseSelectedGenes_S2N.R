# All tumor and control samples
colData_tumor  <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor",])
colData_normal <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal",])
###################################################################################################################################################
# First, I will load the expression table   	
normalized_expression_table<-na.omit(df_reads_count_all_projects[["tpm"]])                                                                       #
                                                                                                                                                 #
# Take the samples ids                                                                                                                           #
sample_stage_I  <-unique(colData_tumor[colData_tumor$stages=="Stage I","sample_id"])                                       #
sample_stage_II <-unique(colData_tumor[colData_tumor$stages=="Stage II","sample_id"])                                      #
sample_stage_III<-unique(colData_tumor[colData_tumor$stages=="Stage III","sample_id"])   
sample_normal   <-unique(colData_normal[,"sample_id"])

# data frame with results
df_results<-data.frame(ENSEMBL=c(), SYMBOL=c(), mean_stage_I=c(), sd_stage_I=c(), log2foldchange_stage_I=c(), pvalue_stage_I=c(), mean_stage_II=c(), sd_stage_II=c(), log2foldchange_stage_II=c(), pvalue_stage_II=c(), mean_stage_III=c(), sd_stage_III=c(), log2foldchange_stage_III=c(), pvalue_stage_III=c())

# List of stage specific genes
stage_specific_genes<-c(unique_stage_I, unique_stage_II, unique_stage_III)
###################################################################################################################################################
# "A total of 4968 up-regulated tumor genes were obtained by comparing all tumor against all normal samples (fdr <=0.05). 
# Among these, 1603 tumor genes are ketpt after filtering for log2foldchange >= 1. Moreover, 6 tumor genes genes whose average expression in normal 
# samples were tpm<=4, because have good signal-to-noise."
# Selected genes

# Vector to store samples labels
df_sample_labels<-data.frame(Samples=unique(c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal)),Tumor=1)

# Storesamples
rownames(df_sample_labels)<-df_sample_labels$Samples
 
# Assert label to samples
df_sample_labels[sample_normal,"Tumor"]<-0


# Sort Table and df_sample_labels
normalized_expression_table<- normalized_expression_table[,c(unique(c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal)))]
df_sample_labels           <- df_sample_labels[c(unique(c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal))),]
#####################################################################################################################################################
ranked_genes<-rank_by_s2n(normalized_expression_table, as.vector(df_sample_labels$Tumor))
#####################################################################################################################################################
