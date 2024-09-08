####################################################################################################################
# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{
	# First, I will load the statistic table   	
	normalized_statistic_table<-read.table(file = paste("/home/felipe/Documents/Cancer_staging/df_statistics_all_projects_",normalization_scheme,".tsv",sep="") , sep = '\t', header = TRUE,fill=TRUE)
	
	# Select only tumor genes
	tumor_genes<-normalized_statistic_table[normalized_statistic_table$tumor_genes == "yes","gene"]
	
	# First, I will load the expression table   	
	normalized_expression_table<-df_reads_count_all_projects[[normalization_scheme]]
	
	# normalized_expression_table
	normalized_expression_table<-normalized_expression_table[tumor_genes,]
	
	# Select only tumor genes
	df_reads_count_all_projects[[normalization_scheme]]<-normalized_expression_table
		
	cat(print(paste("\nNumber of tumor gene after filtering",normalization_scheme,":", length(tumor_genes),"\n")),file=results_files,append=TRUE)
}
