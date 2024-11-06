###################################################################################################################################################
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

# biomarkers
#biomarkers<-data.frame(SYMBOL=c("AKR1B10","GPX2","KRT13","KRT14","KRT16","KRT6B","NTS","S100A7","SPRR1B","SPRR2A"),
#ENSEMBL=c("ENSG00000198074","ENSG00000176153","ENSG00000171401", "ENSG00000186847","ENSG00000186832", "ENSG00000185479", "ENSG00000133636","ENSG00000143556", "ENSG00000169469", "ENSG00000241794"))

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
normalized_expression_table<- na.omit(normalized_expression_table[stage_specific_genes,c(unique(c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal)))])
df_sample_labels           <- df_sample_labels[c(unique(c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal))),]
#####################################################################################################################################################
ranked_genes<-data.frame(rank_by_s2n=rank_by_s2n(normalized_expression_table, as.vector(df_sample_labels$Tumor)))
#####################################################################################################################################################
# Rowmeans for the stage-specific genes
df_rowmeans<-data.frame(RowMeans=(na.omit(rowMeans(normalized_expression_table[,sample_normal]))))

df_rowmeans$Gene<-rownames(df_rowmeans)
ranked_genes$Gene<-rownames(ranked_genes)

# Merge genes
merged_ranked_genes<-merge(df_rowmeans,ranked_genes,by="Gene")

# Set rownames
rownames(merged_ranked_genes)<-merged_ranked_genes$Gene

# Stages collumn
merged_ranked_genes$Stages<-""

merged_ranked_genes[unique_stage_I,"Stages"] <- "Stage I"
merged_ranked_genes[unique_stage_II,"Stages"] <- "Stage II"
merged_ranked_genes[unique_stage_III,"Stages"] <- "Stage III"
#####################################################################################################################################################
# Select genes and collumns
selected_logchange_tumor_control<-list_logchange_tumor_control[["tpm"]][,c("gene", "log2change_all_samples", "pvalue_all_samples", "fdr_all_samples","tumor_genes")]

# Change the colnames
colnames(selected_logchange_tumor_control)<-c("Gene","log2change","pvalue","fdr","tumor")

# Merge ranked genes
merged_ranked_genes_information<-merge(merged_ranked_genes,selected_logchange_tumor_control,by="Gene")

# Re-order genes
merged_ranked_genes_information<-merged_ranked_genes_information[,c("Gene","RowMeans","rank_by_s2n","log2change","pvalue","fdr","tumor","Stages")]
#####################################################################################################################################################
# Path to files of selected_genes                                                                                                             # 
# genes_stages_I
selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

selected_genes_Stage_merged<-rbind(selected_genes_Stage_I_data[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage I"),"Gene"],],
selected_genes_Stage_II_data[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage II"),"Gene"],],
selected_genes_Stage_III_data[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage III"),"Gene"],])

# Set stage
selected_genes_Stage_merged$Stage<-""

selected_genes_Stage_merged[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage I"),"Gene"],"Stage"]<-"Stage I"
selected_genes_Stage_merged[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage II"),"Gene"],"Stage"]<-"Stage II"
selected_genes_Stage_merged[merged_ranked_genes_information[which(merged_ranked_genes_information$Stages=="Stage III"),"Gene"],"Stage"]<-"Stage III"

# Set rownames
rownames(merged_ranked_genes_information)<-merged_ranked_genes_information$Gene
#####################################################################################################################################################
# Set ENSEMBL
df_rowmeans$ENSEMBL <- rownames(df_rowmeans)

# Set biomarkers
merged_ranked_genes_information<-merged_ranked_genes_information[merged_ranked_genes_information$RowMeans <= 4.0,]
selected_genes_Stage_merged<-selected_genes_Stage_merged[rownames(merged_ranked_genes_information),]
#####################################################################################################################################################
# Save TSV file with genes from Stage3
write_tsv(merged_ranked_genes_information, paste(output_dir,"/Statistic_Tumor_Genes.tsv",sep=""))			
write_tsv(selected_genes_Stage_merged, paste(output_dir,"/Statistic_Stage_Specific_Genes.tsv",sep=""))
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


