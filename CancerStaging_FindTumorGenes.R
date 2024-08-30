# A vector with the name of the normalizaton schemes
normalization_schemes <- c("raw","rpkm","fpkm","tpm","tmm")

# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{
  # First, I will load the statistic table   
  normalized_statistic_table<-table<-read.table(file = paste("/home/felipe/Documents/Cancer_staging/df_reads_count_all_projects_",normalization_scheme,".tsv",sep="") , sep = '\t', header = TRUE,fill=TRUE)

  # Second, expression table
  normalized_expression_table<-table<-read.table(file = paste("/home/felipe/Documents/Cancer_staging/df_reads_count_all_projects_",normalization_scheme,".tsv",sep="") , sep = '\t', header = TRUE,fill=TRUE)      
}



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

