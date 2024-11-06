###################################################################################################################################################
# First, I will load the expression table   	
normalized_expression_table<-na.omit(df_reads_count_all_projects[["tpm"]])                                                                       #
                                                                                                                                                 #
# Take the samples ids                                                                                                                           #
sample_stage_I  <-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage I","sample_id"])                                       #
sample_stage_II <-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage II","sample_id"])                                      #
sample_stage_III<-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage III","sample_id"])   
sample_normal   <-unique(merged_data_patient_info[,"sample_id"])

# data frame with results
df_results<-data.frame(ENSEMBL=c(), SYMBOL=c(), mean_stage_I=c(), sd_stage_I=c(), log2foldchange_stage_I=c(), pvalue_stage_I=c(), mean_stage_II=c(), sd_stage_II=c(), log2foldchange_stage_II=c(), pvalue_stage_II=c(), mean_stage_III=c(), sd_stage_III=c(), log2foldchange_stage_III=c(), pvalue_stage_III=c())

# List of stage specific genes
stage_specific_genes<-c(unique_stage_I, unique_stage_II, unique_stage_III)
###################################################################################################################################################
# "A total of 4968 up-regulated tumor genes were obtained by comparing all tumor against all normal samples (fdr <=0.05). 
# Among these, 1603 tumor genes are ketpt after filtering for log2foldchange >= 1. Moreover, 6 tumor genes genes whose average expression in normal 
# samples were tpm<=4, because have good signal-to-noise."
# Selected genes


which(merged_data_patient_info$sample_id %in% colnames(normalized_expression_table))

c(sample_stage_I,sample_stage_II,sample_stage_III)
c(sample_normal)

normalized_expression_table

#####################################################################################################################################################
tpm_stage_I<-selected_genes_Stage_I_data
tpm_stage_II<-selected_genes_Stage_II_data
tpm_stage_III<-selected_genes_Stage_III_data

tpm_stage_I$Stage<-"Stage I"
tpm_stage_II$Stage<-"Stage II"
tpm_stage_III$Stage<-"Stage III"

tpm_stage_all_genes<-rbind(tpm_stage_I,tpm_stage_II,tpm_stage_III)

# Rowmeans for the stage-specific genes
df_rowmeans<-data.frame(RowMeans=(na.omit(rowMeans(normalized_expression_table[rownames(tpm_stage_all_genes),sample_normal]))))
#####################################################################################################################################################
# Set ENSEMBL
df_rowmeans$ENSEMBL <- rownames(df_rowmeans)

# Set biomarkers
biomarkers<-df_rowmeans[df_rowmeans$RowMeans <= 4.0,]
#####################################################################################################################################################
expression_stage_I      <-data.frame(normalized_expression_table[rownames(biomarkers),sample_stage_I])
expression_stage_II     <-data.frame(normalized_expression_table[rownames(biomarkers),sample_stage_II])
expression_stage_III    <-data.frame(normalized_expression_table[rownames(biomarkers),sample_stage_III])
expression_stage_normal <-data.frame(normalized_expression_table[rownames(biomarkers),sample_normal])

expression_stage_I$ENSEMBL<-rownames(expression_stage_I)
expression_stage_II$ENSEMBL<-rownames(expression_stage_II)
expression_stage_III$ENSEMBL<-rownames(expression_stage_III)
expression_stage_normal$ENSEMBL<-rownames(expression_stage_normal)

expression_stage_I          <-melt(expression_stage_I)
expression_stage_II          <-melt(expression_stage_II)
expression_stage_III        <-melt(expression_stage_III)
expression_stage_normal     <-melt(expression_stage_normal)

expression_stage_I$Stages       <-"Tumor"
expression_stage_II$Stages      <-"Tumor"
expression_stage_III$Stages     <-"Tumor"
expression_stage_normal$Stages <-"Control"

expression_all_stages<-rbind(expression_stage_I,expression_stage_II,expression_stage_III,expression_stage_normal)


# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Tumor", "Control"))

# change box plot line colors by groups
p_stage_tumor<-ggplot(expression_all_stages, aes(x=Stages, y=value, fill=Stages)) +   geom_boxplot()+ facet_wrap(~ENSEMBL, ncol = 3, scales="free")+ theme_bw()  

# FindClusters_resolution
png(filename=paste(output_dir,"boplot_selected.png",sep=""), width = 28, height = 28, res=600, units = "cm")
	p_stage_tumor
dev.off()

# Save TSV file with genes from Stage3
write_tsv(na.omit(list_logchange_tumor_control[["tpm"]][rownames(biomarkers),1:4]), paste(output_dir,"/Figure_2_biomarkers_Tumor_Genes.tsv",sep=""))			

################################################################################################################
# Save TSV file with genes from Stage3
write_tsv(tpm_stage_all_genes[rownames(biomarkers),], paste(output_dir,"/Figure_2_biomarkers.tsv",sep=""))			
################################################################################################################
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Stage I", "Control"), c("Stage II", "Control"), c("Stage III", "Control"))

expression_stage_I$Stages       <-"Stage I"
expression_stage_II$Stages      <-"Stage II"
expression_stage_III$Stages     <-"Stage III"
expression_stage_normal$Stages <-"Control"

expression_all_stages<-rbind(expression_stage_I,expression_stage_II,expression_stage_III,expression_stage_normal)

# change box plot line colors by groups
p_stage_stages<-ggplot(expression_all_stages, aes(x=Stages, y=value, fill=Stages)) +   geom_boxplot()+ facet_wrap(~ENSEMBL, ncol = 3, scales="free")+ theme_bw() 

# FindClusters_resolution
png(filename=paste(output_dir,"boplot_selected_per_stage.png",sep=""), width = 32, height = 32, res=600, units = "cm")
	p_stage_stages
dev.off()
