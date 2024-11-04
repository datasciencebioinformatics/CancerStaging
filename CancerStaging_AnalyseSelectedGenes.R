# First, I will load the expression table   	
normalized_expression_table<-na.omit(df_reads_count_all_projects[["tpm"]])

# Take the samples ids
sample_stage_I  <-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage I","sample_id"])                                                          #
sample_stage_II <-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage II","sample_id"])                                                         #
sample_stage_III<-unique(merged_data_patient_info[merged_data_patient_info$stages=="Stage III","sample_id"])   
sample_normal   <-unique(merged_data_patient_info[,"sample_id"])

# biomarkers
#biomarkers<-data.frame(SYMBOL=c("AKR1B10","GPX2","KRT13","KRT14","KRT16","KRT6B","NTS","S100A7","SPRR1B","SPRR2A"),
#ENSEMBL=c("ENSG00000198074","ENSG00000176153","ENSG00000171401", "ENSG00000186847","ENSG00000186832", "ENSG00000185479", "ENSG00000133636","ENSG00000143556", "ENSG00000169469", "ENSG00000241794"))

# data frame with results
df_results<-data.frame(ENSEMBL=c(), SYMBOL=c(), mean_stage_I=c(), sd_stage_I=c(), log2foldchange_stage_I=c(), pvalue_stage_I=c(), mean_stage_II=c(), sd_stage_II=c(), log2foldchange_stage_II=c(), pvalue_stage_II=c(), mean_stage_III=c(), sd_stage_III=c(), log2foldchange_stage_III=c(), pvalue_stage_III=c())

# Selected genes
list_logchange_tumor_control_selected<-list_logchange_tumor_control[["tpm"]][which(list_logchange_tumor_control[["tpm"]]$log2change_all_samples>1),]

# Selected genes
list_logchange_tumor_control_selected<-list_logchange_tumor_control_selected[list_logchange_tumor_control_selected$tumor_genes=="yes",]

# df_rowmeans
df_rowmeans<-data.frame(RowMeans=(na.omit(rowMeans(normalized_expression_table[rownames(list_logchange_tumor_control_selected),sample_normal]))))

# Set ENSEMBL
df_rowmeans$ENSEMBL <- rownames(df_rowmeans)

# Set biomarkers
biomarkers<-df_rowmeans[df_rowmeans$RowMeans <= 4.0,]

# biomarker_ENSEMBL
for (biomarker_ENSEMBL in biomarkers$ENSEMBL)
{
	# SYMBOL
	#gene_symbol<-biomarkers[biomarkers$ENSEMBL==biomarker_ENSEMBL,"SYMBOL"]
	gene_symbol<-"Teste"
	
	# Statistic for stage I
	mean_stage_I           <-mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_I])))
	sd_stage_I             <-sd(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_I]))) 
	log2foldchange_stage_I <-log(mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_I])))/mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal]))),2)
	pvalue_stage_I         <-t.test(x=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_I])), y=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])), paired = FALSE, alternative = "two.sided")$p.value
	
	# Statistic for stage II
	mean_stage_II           <-mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_II])))
	sd_stage_II             <-sd(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_II]))) 
	log2foldchange_stage_II <-log(mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_II])))/mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal]))),2)
	pvalue_stage_II         <-t.test(x=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_II])), y=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])), paired = FALSE, alternative = "two.sided")$p.value
	
	# Statistic for stage III
	mean_stage_III           <-mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_III])))
	sd_stage_III             <-sd(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_III]))) 
	log2foldchange_stage_III <-log(mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_III])))/mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal]))),2)
	pvalue_stage_III         <-t.test(x=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_stage_III])), y=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])), paired = FALSE, alternative = "two.sided")$p.value  

	# Statistic for stage III
	mean_control           <-mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])))
	sd_control             <-sd(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal]))) 
	log2foldchange_control <-log(mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])))/mean(as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal]))),2)
	pvalue_control         <-t.test(x=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])), y=as.vector(t(normalized_expression_table[biomarker_ENSEMBL,sample_normal])), paired = FALSE, alternative = "two.sided")$p.value  
		
	# df_results
	df_results<-rbind(df_results, data.frame(ENSEMBL=biomarker_ENSEMBL, SYMBOL=gene_symbol, mean_stage_I=mean_stage_I, sd_stage_I=sd_stage_I, log2foldchange_stage_I=log2foldchange_stage_I, pvalue_stage_I=pvalue_stage_I, mean_stage_II=mean_stage_II, sd_stage_II=sd_stage_II, log2foldchange_stage_II=log2foldchange_stage_II, pvalue_stage_II=pvalue_stage_II, mean_stage_III=mean_stage_III, sd_stage_III=sd_stage_III, log2foldchange_stage_III=log2foldchange_stage_III, pvalue_stage_III=pvalue_stage_III, mean_control=mean_control, sd_control=sd_control, log2foldchange_control=log2foldchange_control, pvalue_control=pvalue_control))
}
df_results$fdr_stage_I<-p.adjust(df_results$pvalue_stage_I, method="fdr")
df_results$fdr_stage_II<-p.adjust(df_results$pvalue_stage_II, method="fdr")
df_results$fdr_stage_III<-p.adjust(df_results$pvalue_stage_III, method="fdr")

# slected_tumor_genes
slected_tumor_genes<-na.omit(list_logchange_tumor_control[["tpm"]][df_results$ENSEMBL,])

expression_stage_I      <-data.frame(normalized_expression_table[rownames(slected_tumor_genes),sample_stage_I])
expression_stage_II     <-data.frame(normalized_expression_table[rownames(slected_tumor_genes),sample_stage_II])
expression_stage_III    <-data.frame(normalized_expression_table[rownames(slected_tumor_genes),sample_stage_III])
expression_stage_normal <-data.frame(normalized_expression_table[rownames(slected_tumor_genes),sample_normal])

expression_stage_I$ENSEMBL<-rownames(expression_stage_I)
expression_stage_II$ENSEMBL<-rownames(expression_stage_II)
expression_stage_III$ENSEMBL<-rownames(expression_stage_III)
expression_stage_normal$ENSEMBL<-rownames(expression_stage_normal)

expression_stage_I          <-melt(expression_stage_I)
expression_stage_II          <-melt(expression_stage_II)
expression_stage_III        <-melt(expression_stage_III)
expression_stage_normal     <-melt(expression_stage_normal)

expression_stage_I$Stages       <-"Stage I"
expression_stage_II$Stages      <-"Stage II"
expression_stage_III$Stages     <-"Stage III"
expression_stage_normal$Stages <-"Control"


expression_stage_I$Stages       <-"Tumor"
expression_stage_II$Stages      <-"Tumor"
expression_stage_III$Stages     <-"Tumor"
expression_stage_normal$Stages <-"Control"

expression_all_stages<-rbind(expression_stage_I,expression_stage_II,expression_stage_III,expression_stage_normal)

# change box plot line colors by groups
p_stage_tumor<-ggplot(expression_all_stages, aes(x=Stages, y=value, fill=Stages)) +   geom_boxplot()+ facet_wrap(~ENSEMBL, ncol = 3, scales="free")+ theme_bw()

# FindClusters_resolution
png(filename=paste(output_dir,"boplot_selected.png",sep=""), width = 28, height = 14, res=600, units = "cm")
	p_stage_tumor
dev.off()


# Save TSV file with genes from Stage3
write_tsv(na.omit(list_logchange_tumor_control[["tpm"]][df_results$ENSEMBL,]), paste(output_dir,"/Figure_2_biomarkers_Tumor_Genes.tsv",sep=""))			

df_results[df_results$ENSEMBL=="ENSG00000218336",]
list_logchange_tumor_control[["tpm"]]["ENSG00000218336",]


# Save TSV file with genes from Stage3
write_tsv(na.omit(df_results), paste(output_dir,"/Figure_2_biomarkers.tsv",sep=""))			
################################################################################################################
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Stage I", "Control"), c("Stage II", "Control"), c("Stage III", "Control"))
