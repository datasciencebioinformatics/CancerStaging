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
colnames(df_gene_id_symbol)<-c("gene_id","gene_symbol")
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
#############################################################################################################################
# A vector with the name of the normalizaton schemes
normalization_schemes <- c("raw","rpkm","fpkm","tpm","tmm")

# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{
  # First, I will load the statistic table   
  normalized_statistic_table<-table<-read.table(file = paste("/home/felipe/Documents/Cancer_staging/df_statistics_all_projects_",normalization_scheme,".tsv",sep="") , sep = '\t', header = TRUE,fill=TRUE)

  # Second, expression table
  normalized_expression_table<-table<-read.table(file = paste("/home/felipe/Documents/Cancer_staging/df_reads_count_all_projects_",normalization_scheme,".tsv",sep="") , sep = '\t', header = TRUE,fill=TRUE)      

  # Set rownames normalized_statistic_table
  rownames(normalized_statistic_table)<-normalized_statistic_table$gene

  # Set rownames
  rownames(normalized_expression_table)<-normalized_expression_table$gene
  
  # df_gene_id_symbol
  normalized_expression_table<-normalized_expression_table[df_gene_id_symbol$Gene_id,]

  # df_gene_id_symbol
  normalized_statistic_table <-normalized_statistic_table[df_gene_id_symbol$Gene_id,]  

  # Set threshold_normalized
  threshold_normalized <-list_threshold_filters[[normalization_scheme]]

  # Select only the tumor genes
  tumor_genes<-log2change_tumor_control[intersect(which(log2change_tumor_control$fdr_all_samples<=threshold_FDR), which(log2change_tumor_control$log2change_all_samples>=threshold_tumor)),"gene"]  
}
#############################################################################################################################


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
# pvalue_paired           : pvalue_paired tumor/normal paired samples
# fdr_paired               : fdr tumor/normal paired samples


names(list_logchange_tumor_control)

