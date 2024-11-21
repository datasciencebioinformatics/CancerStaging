####################################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm","tpm_calc")
normalization_schemes<-c("tpm")

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{
  #######################################################################################################################################
  # Path to files of selected_genes                                                                                                             # 
  # genes_stages_I
  selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

  rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
  rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
  rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

  selected_genes_Stage_I_gene      <- selected_genes_Stage_I_data$gene
  selected_genes_Stage_II_gene     <- selected_genes_Stage_II_data$gene
  selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene
  #######################################################################################################################################                                                                                                                                     #
  unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
  unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
  unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)
  #############################################################################################################################
  selected_genes_Stage_I_data<-selected_genes_Stage_I_data[unique_stage_I,]
  selected_genes_Stage_II_data<-selected_genes_Stage_II_data[unique_stage_II,]
  selected_genes_Stage_III_data<-selected_genes_Stage_III_data[unique_stage_III,]

  colnames(selected_genes_Stage_I_data)[1]<-c("ensembl_gene_id")
  colnames(selected_genes_Stage_II_data)[1]<-c("ensembl_gene_id")
  colnames(selected_genes_Stage_III_data)[1]<-c("ensembl_gene_id")

  genes_rankData_stage_I     <- getBM(filters= "ensembl_genwrite_tsv(na.omit(genes_rankData_stage_all_genes[selected_genes_Stage_I_gene,]),   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_I",".tsv",sep=""))e_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_I_data),mart=mart)
  genes_rankData_stage_II    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_II_data),mart=mart)
  genes_rankData_stage_III   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_III_data),mart=mart)

  genes_rankData_stage_I<-merge(selected_genes_Stage_I_data,genes_rankData_stage_I,by="ensembl_gene_id")
  genes_rankData_stage_II<-merge(selected_genes_Stage_II_data,genes_rankData_stage_II,by="ensembl_gene_id")
  genes_rankData_stage_III<-merge(selected_genes_Stage_III_data,genes_rankData_stage_III,by="ensembl_gene_id")

  write_tsv(na.omit(genes_rankData_stage_all_genes[selected_genes_Stage_I_gene,]),   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_I",".tsv",sep=""))
  write_tsv(genes_rankData_stage_I,     paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_I",".tsv",sep=""))
  write_tsv(genes_rankData_stage_II,   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_II",".tsv",sep=""))
  write_tsv(genes_rankData_stage_III,  paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_III",".tsv",sep=""))  

    
}










# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{
  #######################################################################################################################################
  # Path to files of selected_genes                                                                                                             # 
  # genes_stages_I
  selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

  rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
  rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
  rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

  selected_genes_Stage_I_gene      <- selected_genes_Stage_I_data$gene
  selected_genes_Stage_II_gene     <- selected_genes_Stage_II_data$gene
  selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene
  #######################################################################################################################################                                                                                                                                     #
  unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
  unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
  unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)
  #############################################################################################################################
  selected_genes_Stage_I_data<-selected_genes_Stage_I_data[unique_stage_I,]
  selected_genes_Stage_II_data<-selected_genes_Stage_II_data[unique_stage_II,]
  selected_genes_Stage_III_data<-selected_genes_Stage_III_data[unique_stage_III,]

  colnames(selected_genes_Stage_I_data)[1]<-c("ensembl_gene_id")
  colnames(selected_genes_Stage_II_data)[1]<-c("ensembl_gene_id")
  colnames(selected_genes_Stage_III_data)[1]<-c("ensembl_gene_id")

  genes_rankData_stage_I     <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_I_data),mart=mart)
  genes_rankData_stage_II    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_II_data),mart=mart)
  genes_rankData_stage_III   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_III_data),mart=mart)

  genes_rankData_stage_I<-merge(selected_genes_Stage_I_data,genes_rankData_stage_I,by="ensembl_gene_id")
  genes_rankData_stage_II<-merge(selected_genes_Stage_II_data,genes_rankData_stage_II,by="ensembl_gene_id")
  genes_rankData_stage_III<-merge(selected_genes_Stage_III_data,genes_rankData_stage_III,by="ensembl_gene_id")

  write_tsv(genes_rankData_stage_I,     paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_I",".tsv",sep=""))
  write_tsv(genes_rankData_stage_II,   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_II",".tsv",sep=""))
  write_tsv(genes_rankData_stage_III,  paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","unique_stage_III",".tsv",sep=""))  

    
}
#######################################################################################################################################
normalization_scheme<-"tpm"
#######################################################################################################################################
# Path to files of selected_genes                                                                                                             # 
# genes_stages_I
selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

selected_genes_Stage_I_gene      <- selected_genes_Stage_I_data$gene
selected_genes_Stage_II_gene     <- selected_genes_Stage_II_data$gene
selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene
#######################################################################################################################################                                                                                                                                     #
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
write_tsv(genes_rankData_stage_all_genes_merged,   paste(output_dir,"/tumor_genes_statistics",normalization_scheme,".tsv",sep=""))

Table3<-genes_rankData_stage_all_genes[genes_rankData_stage_all_genes$FC>=50 & genes_rankData_stage_all_genes$Mean_normal<=10,]

# Save TSV file with genes from Stage3
write_tsv(Table3, paste(output_dir,"/Table3.tsv",sep=""))



# Check tumor table
table_tumor_genes  <-na.omit(genes_rankData_stage_all_genes)
table_stage_I_genes<-na.omit(genes_rankData_stage_all_genes[selected_genes_Stage_I_gene,])
table_stage_II_genes<-na.omit(genes_rankData_stage_all_genes[selected_genes_Stage_II_gene,])
table_stage_III_genes<-na.omit(genes_rankData_stage_all_genes[selected_genes_Stage_III_gene,])

table_stage_I_genes$FC<-table_stage_I_genes$FC_Stage_I
table_stage_II_genes$FC<-table_stage_II_genes$FC_Stage_II
table_stage_III_genes$FC<-table_stage_III_genes$FC_Stage_III

table_stage_I_genes<-table_stage_I_genes[,c("ensembl_gene_id","Mean_normal","sd_normal","Mean_tumor","sd_tumor","FC","Mean_Stage_I","sd_Stage_I","Mean_Stage_II","sd_Stage_II","Mean_Stage_III","sd_Stage_III","entrezgene_accession","entrezgene_id","hgnc_symbol","description")]
table_stage_II_genes<-table_stage_II_genes[,c("ensembl_gene_id","Mean_normal","sd_normal","Mean_tumor","sd_tumor","FC","Mean_Stage_I","sd_Stage_I","Mean_Stage_II","sd_Stage_II","Mean_Stage_III","sd_Stage_III","entrezgene_accession","entrezgene_id","hgnc_symbol","description")]
table_stage_III_genes<-table_stage_III_genes[,c("ensembl_gene_id","Mean_normal","sd_normal","Mean_tumor","sd_tumor","FC","Mean_Stage_I","sd_Stage_I","Mean_Stage_II","sd_Stage_II","Mean_Stage_III","sd_Stage_III","entrezgene_accession","entrezgene_id","hgnc_symbol","description")]

write_tsv(table_tumor_genes,   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_tumor_genes",".tsv",sep=""))
write_tsv(na.omit(table_stage_I_genes),   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_I",".tsv",sep=""))
write_tsv(na.omit(table_stage_II_genes),  paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_II",".tsv",sep=""))
write_tsv(na.omit(table_stage_III_genes), paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_III",".tsv",sep=""))
