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

# Take the gene os each stage
stage_I_genes  <-list_stage_specific_genes[["tpm_stage_I"]][list_stage_specific_genes[["tpm_stage_I"]]$Category=="Per stage genes",]
stage_II_genes  <-list_stage_specific_genes[["tpm_stage_II"]][list_stage_specific_genes[["tpm_stage_II"]]$Category=="Per stage genes",]
stage_III_genes  <-list_stage_specific_genes[["tpm_stage_III"]][list_stage_specific_genes[["tpm_stage_III"]]$Category=="Per stage genes",]

rownames(stage_I_genes)<-stage_I_genes$gene
rownames(stage_II_genes)<-stage_II_genes$gene
rownames(stage_III_genes)<-stage_III_genes$gene

colnames(stage_I_genes)[1]<-c("ensembl_gene_id")
colnames(stage_II_genes)[1]<-c("ensembl_gene_id")
colnames(stage_III_genes)[1]<-c("ensembl_gene_id")

genes_rankData_stage_I     <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_I_data),mart=mart)
genes_rankData_stage_II    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_II_data),mart=mart)
genes_rankData_stage_III   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol","description"),values=rownames(selected_genes_Stage_III_data),mart=mart)

genes_rankData_stage_I<-merge(stage_I_genes,genes_rankData_stage_I,by="ensembl_gene_id")
genes_rankData_stage_II<-merge(stage_II_genes,genes_rankData_stage_II,by="ensembl_gene_id")
genes_rankData_stage_III<-merge(stage_III_genes,genes_rankData_stage_III,by="ensembl_gene_id")

genes_rankData_stage_I <- genes_rankData_stage_I[match(unique(genes_rankData_stage_I$ensembl_gene_id), genes_rankData_stage_I$ensembl_gene_id),]
genes_rankData_stage_II <- genes_rankData_stage_I[match(unique(genes_rankData_stage_II$ensembl_gene_id), genes_rankData_stage_II$ensembl_gene_id),]
genes_rankData_stage_III <- genes_rankData_stage_I[match(unique(genes_rankData_stage_III$ensembl_gene_id), genes_rankData_stage_III$ensembl_gene_id),]

write_tsv(genes_rankData_stage_I,     paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_I",".tsv",sep=""))
write_tsv(genes_rankData_stage_II,   paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_II",".tsv",sep=""))
write_tsv(genes_rankData_stage_III,  paste(output_dir,"/FindStageSpecificGenes_annotated",normalization_scheme,"_","all_stage_III",".tsv",sep=""))  

