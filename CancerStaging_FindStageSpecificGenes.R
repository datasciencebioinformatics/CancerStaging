#######################################################################################################################################
# All tumor and control samples
colData_tumor  <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor",])
colData_normal <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal",])
                                                                                                                                     #
# Vector with each stage                                                                                                              #
stages_str<-c("stage_I","stage_II","stage_III")                                                                                       #
                                                                                                                                      #
# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-unique(colData_tumor[colData_tumor$stages=="Stage I","sample_id"])                                                          #
sample_stage_II <-unique(colData_tumor[colData_tumor$stages=="Stage II","sample_id"])                                                         #
sample_stage_III<-unique(colData_tumor[colData_tumor$stages=="Stage III","sample_id"])                                                        #
sample_normal   <-unique(colData_normal[,"sample_id"])                                                                                        #
#######################################################################################################################################
df_table_comparisson=rbind(data.frame(Stage_i="sample_stage_I",Stage_ii="sample_normal"),                                             #
data.frame(Stage_i="sample_stage_II",Stage_ii="sample_normal"),                                                                       # 
data.frame(Stage_i="sample_stage_III",Stage_ii="sample_normal"))                                                                      #
#######################################################################################################################################
list_of_comparisson=list(sample_stage_I=sample_stage_I,sample_stage_II=sample_stage_II,sample_stage_III=sample_stage_III, sample_normal=sample_normal)
#######################################################################################################################################
list_stage_specific_genes<-c()
####################################################################################################################
# for each  normalization scheme
for (normalization_scheme in names(df_reads_count_all_projects))
{
	# First, I will load the statistic table   	
	normalized_statistic_table<-list_logchange_tumor_control[[normalization_scheme]]
	
	# Select only tumor genes
	normalized_statistic_table<-normalized_statistic_table[normalized_statistic_table$tumor_genes == "yes",]
	
	# First, I will load the expression table   	
	normalized_expression_table<-na.omit(df_reads_count_all_projects[[normalization_scheme]])
	
	# for each pair of stage.
	for (comparisson_index in rownames(df_table_comparisson))
	{	
		# Stages
		Stage_i          <-df_table_comparisson[comparisson_index,"Stage_i"]
		Stage_ii         <-df_table_comparisson[comparisson_index,"Stage_ii"]
		
		# Take gens of corresponding stage
		DE_genes        <- unique(normalized_statistic_table$gene)
		
		# Take samples of each stage
		Stage_i_samples         =list_of_comparisson[[Stage_i]]
		Stage_ii_samples        =list_of_comparisson[[Stage_ii]]	
		
		# Take RPKM of genes from samples of each stage
		Stages_i_samples_expr         <-na.omit(normalized_expression_table[DE_genes,Stage_i_samples])
		Stages_ii_samples_expr        <- na.omit(normalized_expression_table[DE_genes,Stage_ii_samples])
		
		####################################################################################################################
		# folchange=Expr(Stage i)/Expr(Stage ii and II)
		#folchange=rowMeans(Stage_i_samples_expr)/rowMeans(Stages_ii_sample_expr)	
		#log2change=log(folchange,2)	
		
		# Log2foldchange
		LOG_CONSTANT=0.001
		log2change=log( (rowMeans(Stages_i_samples_expr+LOG_CONSTANT)/rowMeans(Stages_ii_samples_expr+LOG_CONSTANT)),2)	
		
		# log2change data
		log2change_Stage_i=na.omit(data.frame(gene=names(log2change),log2change=log2change))
		
		# First, the log2foldchane tumor/normal samples is used
		log2change_Stage_i$Category<-"insignificant"
		
		# First, the log2foldchane tumor/normal samples is used
		log2change_Stage_i$Pvalue<-1
		
		# For each genes in the tabe
		for (gene in log2change_Stage_i$gene)
		{
			# Take expression only for the gene
			Stage_i_gene_expr<-Stages_i_samples_expr[gene,]
			Stage_ii_gene_expr<-Stages_ii_samples_expr[gene,]
			
			# Filter by threshold_filters
			Stage_i_gene_expr<-Stage_i_gene_expr[Stage_i_gene_expr > list_threshold_filters[[normalized_table_names]]]
		
			# Filter by threshold_filters
			Stage_ii_gene_expr<-Stage_ii_gene_expr[Stage_ii_gene_expr > list_threshold_filters[[normalized_table_names]]]	  
					
			# Take p-value				
			out <- tryCatch(log2change_Stage_i[gene,"Pvalue"]<-t.test(x=as.numeric(Stage_i_gene_expr), y=as.numeric(Stage_ii_gene_expr), paired = FALSE, alternative = "two.sided")$p.value, error = function(e) NULL)
		}
		# FRD 
		log2change_Stage_i$FDR<-p.adjust(log2change_Stage_i$Pvalue, method="fdr")
		
		# Categorize genes if log2foldchange >= threshold_
		log2change_Stage_i[intersect(which(log2change_Stage_i$FDR<=threshold_FDR), which(log2change_Stage_i$log2change>=threshold_stage)),"Category"]<-paste("Per stage genes", sep="")
		#log2change_Stage_i[which(log2change_Stage_i$FDR<=0.05),"Category"]<-paste("Per stage genes", sep="")
		
		# Selected genes based on FDR, log2foldchange
		selected_genes<-log2change_Stage_i[intersect(which(log2change_Stage_i$FDR<=threshold_FDR), which(log2change_Stage_i$log2change>=threshold_stage)),"gene"]	
		#selected_genes<-log2change_Stage_i[ which(log2change_Stage_i$log2change>=threshold_stage),"gene"]	
		####################################################################################################################	
		# Save TSV file with genes from Stage3
		write_tsv(na.omit(log2change_Stage_i[selected_genes,]), paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_",Stage_i,".tsv",sep=""))			
		####################################################################################################################		
		list_stage_specific_genes[[paste(normalization_scheme,"_",substring(Stage_i,8,20),sep="")]]<-log2change_Stage_i
		####################################################################################################################	
		cat(print(paste("\nStage-specific genes ",Stage_i, " : " ,normalization_scheme," : ",length(selected_genes),"\n",sep="")),file=results_files,append=TRUE)
	}			
	# Save TSV file with genes from Stage3
	write_tsv(normalized_statistic_table, paste(output_dir,"df_statistics_all_projects_",normalization_scheme,".tsv",sep=""))	
	
}
save(list_stage_specific_genes, file = paste(output_dir,"/","StageSpecificGenes.RData",sep=""))
print("\nCancerStaging_FindStageSpecificGenes")
