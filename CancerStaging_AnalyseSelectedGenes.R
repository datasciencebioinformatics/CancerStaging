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
biomarkers<-

data.frame(SYMBOL=c("AKR1B10","GPX2","KRT13","KRT14","KRT16","KRT6B","NTF","S100A7","SPRR1B","SPRR2A"),
ENSEMBL=c("ENSG00000198074","ENSG00000176153","ENSG00000171401", "ENSG00000186847","ENSG00000186832", "ENSG00000185479", "ENSG00000267019","ENSG00000143556", "ENSG00000169469", "ENSG00000241794"))

for (biomarker in biomarkers)
{
  # convenrt 
  EnsemblToUniprotKBconversionList_data
  
  normalized_expression_table[,sample_stage_I]
  
}

