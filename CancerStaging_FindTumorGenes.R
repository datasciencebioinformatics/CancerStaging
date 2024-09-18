############################################################################################################################3
gene_ids_file                    <- "/home/felipe/Documents/Cancer_staging/tables/gene_name.txt"
gene_name_file                   <- "/home/felipe/Documents/Cancer_staging/tables/gene_ids.txt"
EnsemblToUniprotKBconversionList <- "/home/felipe/Documents/Cancer_staging/EnsemblToUniprotKBconversionList.txt"
Table_1                          <- "/home/felipe/Documents/Cancer_staging/Table_1.tsv"
ListGenesInteratoma              <- "/home/felipe/Documents/Cancer_staging/ListGenesInteratoma.tsv"
############################################################################################################################3
# Load gene data (name and id)
gene_ids_data        <-read.table(file = gene_ids_file, sep = '\t', header = FALSE,fill=TRUE)      
gene_name_data       <-read.table(file = gene_name_file, sep = '\t', header = FALSE,fill=TRUE)      
Table1_data          <-read.table(file = Table_1, sep = '\t', header = FALSE,fill=TRUE)    
Table2_interactoma   <-read.table(file = ListGenesInteratoma, sep = '\t', header = TRUE,fill=TRUE)    
############################################################################################################################
# Colnames gene_id and gene_symbol
colnames(Table1_data)<-c("gene_id","gene_symbol")

# Put both in a data
df_gene_id_symbol<-data.frame(gene_id=gene_ids_data,gene_symbol=gene_name_data)

# Rename collumns
colnames(df_gene_id_symbol)<-c("gene_symbol","gene_id")
#############################################################################################################################
# Set Gene_id
df_gene_id_symbol$Gene_id <- ""

# For each gene
for (gene in df_gene_id_symbol$gene_id)
{
  # Take the gene id
  gene_id<-substring(gene,1,last=15)

  # df_gene_id_symbol
  df_gene_id_symbol[which(grepl( gene_id,df_gene_id_symbol$gene_id, fixed = TRUE)),"Gene_id"]<-gene_id
}

# Merge Table2_interactoma and df_gene_id_symbol
merge_interactome_gene_symbol<-merge(x=Table2_interactoma, y=df_gene_id_symbol, by = "gene_symbol")
#####################################################################################################################
# Find tumor genes by padj and log2foldchange
# In this table, there are the statistics for each of the normalization scheme.
# The statistics compare the tumor against normal samples in two way,
# first all turmor against all normals, second, only the paired tumor against the paired normal.
# The name of the collumbns are:
# gene                     : ENSEMBL
# log2change_all_samples   : log2change tumor/normal all tumor samples/all control samples
# pvalue_all_samples       : pvalue tumor/normal all tumor samples/all control samples
# fdr_all_samples          : fdr tumor/normal all tumor samples/all control samples
# log2change_paired        : log2change_paired tumor/normal paired samples
# pvalue_paired            : pvalue_paired tumor/normal paired samples
# fdr_paired               : fdr tumor/normal paired samples
#############################################################################################################################
# for each  normalization scheme
for (normalization_scheme in names(df_reads_count_all_projects))
{	
	# First, I will load the statistic table   	
	normalized_statistic_table<-list_logchange_tumor_control[[normalization_scheme]]
		
	# Set rownames normalized_statistic_table
	rownames(normalized_statistic_table)<-normalized_statistic_table$gene
	
	# df_gene_id_symbol
	normalized_statistic_table <-normalized_statistic_table[unique(merge_interactome_gene_symbol$Gene_id),]
	
	# Set threshold_normalized
	threshold_normalized <-list_threshold_filters[[normalization_scheme]]
	
	# Set tumor genes field
	normalized_statistic_table$tumor_genes <- "no"
	
	# Select only the tumor genes
	normalized_statistic_table[intersect(which(normalized_statistic_table$fdr_all_samples<=threshold_FDR), which(normalized_statistic_table$log2change_all_samples>=threshold_tumor)),"tumor_genes"]  <- "yes"

	# First, I will load the statistic table   	
	normalized_statistic_table<-na.omit(normalized_statistic_table)

	# First, I will load the statistic table   	
	list_logchange_tumor_control[[normalization_scheme]]<-normalized_statistic_table	

	cat(print(paste("\nNumber of tumor gene :", paste(normalization_scheme," : ",sum(normalized_statistic_table$tumor_genes=="yes"),"",sep=" "))),file=results_files,append=TRUE)
		
	# Save TSV file with genes from Stage3
	write_tsv(normalized_statistic_table, paste(output_dir,"df_statistics_all_projects_",normalization_scheme,".tsv",sep=""))			  
}
print("\nCancerStaging_FindTumorGenes")
