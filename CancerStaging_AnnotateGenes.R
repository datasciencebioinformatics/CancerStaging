########################################################################################################################
# Set rowanames
rownames(df_FC)<-rownames(df_FC)

# Set colnames
colnames(df_FC)[1]<-c("ensembl_gene_id")

# df_FC
genes_df_FC    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(df_FC),mart=mart)

#  Merge data.frame
genes_rankData_stage_all_genes<-merge(df_FC,genes_df_FC,by="ensembl_gene_id")

# Take only the first occurance of a  
genes_rankData_stage_all_genes <- genes_rankData_stage_all_genes[match(unique(genes_rankData_stage_all_genes$ensembl_gene_id), genes_rankData_stage_all_genes$ensembl_gene_id),]

# Set rownames
rownames(genes_rankData_stage_all_genes)<-genes_rankData_stage_all_genes$ensembl_gene_id

# Set colnames
colnames(Interactomes_GC3_T2_merged)[4]<-"ensembl_gene_id"

# Set genes_rankData_stage_all_genes
genes_rankData_stage_all_genes_merged<-merge(Interactomes_GC3_T2_merged, genes_rankData_stage_all_genes,by="ensembl_gene_id")

# genes_rankData_stage_all_genes_merged
write_tsv(genes_rankData_stage_all_genes_merged,   paste(output_dir,"/Table_S3.tsv",sep=""))

Table3<-genes_rankData_stage_all_genes[genes_rankData_stage_all_genes$FC>=50 & genes_rankData_stage_all_genes$Mean_normal<=10,]

# Save TSV file with genes from Stage3
write_tsv(Table3, paste(output_dir,"/Table8.tsv",sep=""))
