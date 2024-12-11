###########################################################################################################
# Restore the object                                                                                      #
Interactomes_GC3_T2_merged <-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))    #
normalization_schemes      <-readRDS(file = paste(output_dir,"normalization_schemes.rds",sep=""))         #
df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))   #
list_of_comparisson        <-readRDS(file = paste(output_dir,"list_of_comparisson.rds",sep=""))           #
###########################################################################################################
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

# Format table 8
Table8<-genes_rankData_stage_all_genes_merged[genes_rankData_stage_all_genes_merged$FC>=50 & genes_rankData_stage_all_genes_merged$Mean_normal<=10,c("hgnc_symbol","Mean_normal","sd_normal","Mean_tumor","sd_tumor","FC","Conections")]

# Save TSV file with genes from Stage3
write_tsv(Table8, paste(output_dir,"/Table8.tsv",sep=""))
#########################################################################################################
# Format table S3
genes_rankData_stage_all_genes_merged$FC5 <-"NO"
genes_rankData_stage_all_genes_merged$FC10<-"NO"
genes_rankData_stage_all_genes_merged$FC20<-"NO"
genes_rankData_stage_all_genes_merged$FC30<-"NO"
genes_rankData_stage_all_genes_merged$FC40<-"NO"
genes_rankData_stage_all_genes_merged$FC50<-"NO"

genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=5),"FC5"]   <-"YES"
genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=10),"FC10"] <-"YES"
genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=20),"FC20"] <-"YES"
genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=30),"FC30"] <-"YES"
genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=40),"FC40"] <-"YES"
genes_rankData_stage_all_genes_merged[which(genes_rankData_stage_all_genes_merged$FC>=50),"FC50"] <-"YES"
#########################################################################################################
# genes_rankData_stage_all_genes_merged
write_tsv(genes_rankData_stage_all_genes_merged,   paste(output_dir,"/Table_S3.tsv",sep=""))
#########################################################################################################

