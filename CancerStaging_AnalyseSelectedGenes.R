# First, I will load the expression table   	
normalized_expression_table<-na.omit(df_reads_count_all_projects[["tpm"]])

# Take the samples ids
sample_stage_I  <-unique(colData_tumor[colData_tumor$stages=="Stage I","sample_id"])                                                          #
sample_stage_II <-unique(colData_tumor[colData_tumor$stages=="Stage II","sample_id"])                                                         #
sample_stage_III<-unique(colData_tumor[colData_tumor$stages=="Stage III","sample_id"])   
sample_normal   <-unique(colData_normal[,"sample_id"])

LOG_CONSTANT=0.001
log2change=log( (rowMeans(Stages_i_samples_expr+LOG_CONSTANT)/rowMeans(Stages_ii_samples_expr+LOG_CONSTANT)),2)	

# biomarkers
biomarkers<-data.frame(SYMBOL=c("AKR1B10","GPX2","KRT13","KRT14","KRT16","KRT6B","NTS","S100A7","SPRR1B","SPRR2A"),
ENSEMBL=c("ENSG00000198074","ENSG00000176153","ENSG00000171401", "ENSG00000186847","ENSG00000186832", "ENSG00000185479", "ENSG00000267019","ENSG00000143556", "ENSG00000169469", "ENSG00000241794"))

# data frame with results
df_results<-data.frame(ENSEMBL=c(), SYMBOL=c(), mean_stage_I=c(), sd_stage_I=c(), log2foldchange_stage_I=c(), pvalue_stage_I=c(), mean_stage_II=c(), sd_stage_II=c(), log2foldchange_stage_II=c(), pvalue_stage_II=c(), mean_stage_III=c(), sd_stage_III=c(), log2foldchange_stage_III=c(), pvalue_stage_III=c())

for (biomarker_ENSEMBL in biomarkers$ENSEMBL)
{
	# SYMBOL
	gene_symbol<-biomarkers[biomarkers$ENSEMBL==biomarker_ENSEMBL,"SYMBOL"]
	
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
	
	# df_results
	df_results<-rbind(df_results, data.frame(ENSEMBL=biomarker_ENSEMBL, SYMBOL=gene_symbol, mean_stage_I=mean_stage_I, sd_stage_I=sd_stage_I, log2foldchange_stage_I=log2foldchange_stage_I, pvalue_stage_I=pvalue_stage_I, mean_stage_II=mean_stage_II, sd_stage_II=sd_stage_II, log2foldchange_stage_II=log2foldchange_stage_II, pvalue_stage_II=pvalue_stage_II, mean_stage_III=mean_stage_III, sd_stage_III=sd_stage_III, log2foldchange_stage_III=log2foldchange_stage_III, pvalue_stage_III=pvalue_stage_III))
}


