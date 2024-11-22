library("msigdb")
library("msigdbr")
library("fgsea")

# Load mart tables
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Install msigdb packages.
# use the custom accessor to select a specific version of MSigDB
msigdb.hs = getMsigdb(org = 'hs', id = 'EZID', version = '7.4')

# Use also the msigdb
pathwaysDF <- msigdbr("human", category="H")

# Retrieve the pathways
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)



# retrieeve the hallmarks gene sets
# 50 hallmarks
hallmarks_gene_set<-subsetCollection(msigdb.hs, 'h')

# First, I will load the expression table   	
# Take the expression of genes from sameples of each stage
#expr_stage_I<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_I_gene,sample_stage_I])
#expr_stage_II<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_II_gene,sample_stage_II])
#expr_stage_III<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_III_gene,sample_stage_III])

# Omit lines with NA
expr_stage_I<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_I,sample_stage_I])
expr_stage_II<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_II,sample_stage_II])
expr_stage_III<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_III,sample_stage_III])

# Take for each ensembl_gene_id the entrezgene_accession, entrezgene_id, hgnc_symbol
genes_rankData_stage_I     <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_I),mart=mart)
genes_rankData_stage_II    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_II),mart=mart)
genes_rankData_stage_III   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_III),mart=mart)

# Next, Take the annotation for each genes
genes_rankData_stage_I   <-genes_rankData_stage_I[genes_rankData_stage_I$ensembl_gene_id %in% rownames(expr_stage_I),]
genes_rankData_stage_II  <-genes_rankData_stage_II[genes_rankData_stage_II$ensembl_gene_id %in% rownames(expr_stage_II),]
genes_rankData_stage_III <-genes_rankData_stage_III[genes_rankData_stage_III$ensembl_gene_id %in% rownames(expr_stage_III),]

# Take the first occurance of each ensembl_gene_id 
genes_rankData_stage_I   <- genes_rankData_stage_I[match(unique(genes_rankData_stage_I$ensembl_gene_id), genes_rankData_stage_I$ensembl_gene_id),]
genes_rankData_stage_II  <- genes_rankData_stage_II[match(unique(genes_rankData_stage_II$ensembl_gene_id), genes_rankData_stage_II$ensembl_gene_id),]
genes_rankData_stage_III <- genes_rankData_stage_III[match(unique(genes_rankData_stage_III$ensembl_gene_id), genes_rankData_stage_III$ensembl_gene_id),]

# Set the rownames as the ensembl_gene_id
rownames(genes_rankData_stage_I)<-genes_rankData_stage_I$ensembl_gene_id
rownames(genes_rankData_stage_II)<-genes_rankData_stage_II$ensembl_gene_id
rownames(genes_rankData_stage_III)<-genes_rankData_stage_III$ensembl_gene_id

# Set the rownames entrezgene_id
rownames(expr_stage_I)   <- genes_rankData_stage_I[rownames(expr_stage_I),"entrezgene_id"]
rownames(expr_stage_II)  <- genes_rankData_stage_II[rownames(expr_stage_II),"entrezgene_id"]
rownames(expr_stage_III) <- genes_rankData_stage_III[rownames(expr_stage_III),"entrezgene_id"]

# Compute the geseca
geseca_Stage_I   <- data.frame(geseca(pathways, expr_stage_I))
geseca_Stage_II  <- data.frame(geseca(pathways, expr_stage_II))
geseca_Stage_III <- data.frame(geseca(pathways, expr_stage_III))

# Filter padj
geseca_Stage_I  <-geseca_Stage_I[geseca_Stage_I$padj<=0.05,]
geseca_Stage_II <-geseca_Stage_II[geseca_Stage_II$padj<=0.05,]
geseca_Stage_III<-geseca_Stage_III[geseca_Stage_III$padj<=0.05,]

geseca_Stage_I$Stage<-"Stage I"
geseca_Stage_II$Stage<-"Stage II"
geseca_Stage_III$Stage<-"Stage III"
#############################################################################################################################################################################################3

# Combine the queries for the three stages
#merged_pathways<-merge(merge(geseca_Stage_I,geseca_Stage_II,by="pathway", all=TRUE),geseca_Stage_III,by="pathway", all=TRUE)
merge(geseca_Stage_II,geseca_Stage_III,by="pathway", all=TRUE)


# Table with results
results_hallmark<-data.frame(hallmarks=hallmarks,
             genes_n_Stage_I=genes_n_Stage_I,
             genes_per_Stage_I=genes_percentage_Stage_I,
             genes_padj_Stage_I=genes_padj_Stage_I,
             genes_n_Stage_II=genes_n_Stage_II,
             genes_per_Stage_II=genes_percentage_Stage_II,
             genes_padj_Stage_II=genes_padj_Stage_II,
             genes_n_Stage_III=genes_n_Stage_III,
             genes_per_Stage_III=genes_percentage_Stage_III,
             genes_padj_Stage_III=genes_padj_Stage_III,                                                        
             symbol_Stage_I=symbol_Stage_I,
             symbol_Stage_II=symbol_Stage_II,
             symbol_Stage_III=symbol_Stage_III)

# For each hallmark
for (hallmarks in names(pathways))
{
  # Take each genes
  symbol_Stage_I<-paste(genes_rankData_stage_I[rownames(expr_stage_I) %in% pathways[[hallmarks]],"hgnc_symbol"],collapse=" , ")
  symbol_Stage_II<-paste(genes_rankData_stage_II[rownames(expr_stage_II) %in% pathways[[hallmarks]],"hgnc_symbol"],collapse=" , ")
  symbol_Stage_III<-paste(genes_rankData_stage_III[rownames(expr_stage_III) %in% pathways[[hallmarks]],"hgnc_symbol"],collapse=" , ")

  # Take number of genes from this ptahway on stage I
  genes_Stage_I<-paste(rownames(expr_stage_I)[rownames(expr_stage_I) %in% pathways[[hallmarks]]],collapse=" , ")
  genes_Stage_II<-paste(rownames(expr_stage_II)[rownames(expr_stage_II) %in% pathways[[hallmarks]]],collapse=" , ")
  genes_Stage_III<-paste(rownames(expr_stage_III)[rownames(expr_stage_III) %in% pathways[[hallmarks]]],collapse=" , ")

  # Take number of genes from this ptahway on stage I
  genes_n_Stage_I<-sum(rownames(expr_stage_I) %in% pathways[[hallmarks]])
  genes_n_Stage_II<-sum(rownames(expr_stage_II) %in% pathways[[hallmarks]])
  genes_n_Stage_III<-sum(rownames(expr_stage_III) %in% pathways[[hallmarks]])

 # Take number of genes from this ptahway on stage I
  genes_percentage_Stage_I<-genes_n_Stage_I/length(rownames(expr_stage_I))
  genes_percentage_Stage_II<-genes_n_Stage_II/length(rownames(expr_stage_II))
  genes_percentage_Stage_III<-genes_n_Stage_III/length(rownames(expr_stage_III))

  # Take number of genes from this ptahway on stage I
  padj_Stage_I<-genes_n_Stage_I/length(rownames(expr_stage_I))
  padj_Stage_II<-genes_n_Stage_II/length(rownames(expr_stage_II))
  padj_Stage_III<-genes_n_Stage_III/length(rownames(expr_stage_III))

  # Take the pad of 
  genes_padj_Stage_I   <-geseca_Stage_I[geseca_Stage_I$pathway==hallmarks,"padj"]
  genes_padj_Stage_II  <-geseca_Stage_II[geseca_Stage_II$pathway==hallmarks,"padj"]
  genes_padj_Stage_III <-geseca_Stage_III[geseca_Stage_III$pathway==hallmarks,"padj"]

  if(length(genes_padj_Stage_I)==0)
    genes_padj_Stage_I<-1
  if(length(genes_padj_Stage_II)==0)
    genes_padj_Stage_II<-1
  if(length(genes_padj_Stage_III)==0)
    genes_padj_Stage_III<-1  
  
  # If pathway has representation
  if(genes_padj_Stage_I+genes_padj_Stage_II+genes_padj_Stage_III<3)
  {
    # Table with results
    results_hallmark<-rbind(results_hallmark,data.frame(hallmarks=hallmarks,
             genes_n_Stage_I=genes_n_Stage_I,
             genes_per_Stage_I=genes_percentage_Stage_I,
             genes_padj_Stage_I=genes_padj_Stage_I,
             genes_n_Stage_II=genes_n_Stage_II,
             genes_per_Stage_II=genes_percentage_Stage_II,
             genes_padj_Stage_II=genes_padj_Stage_II,
             genes_n_Stage_III=genes_n_Stage_III,
             genes_per_Stage_III=genes_percentage_Stage_III,
             genes_padj_Stage_III=genes_padj_Stage_III,                                                        
             symbol_Stage_I=symbol_Stage_I,
             symbol_Stage_II=symbol_Stage_II,
             symbol_Stage_III=symbol_Stage_III))    
  }  
}
# results_hallmark
results_hallmark<-results_hallmark[,c("hallmarks","symbol_Stage_I","genes_n_Stage_I","genes_per_Stage_I","genes_padj_Stage_I","symbol_Stage_II","genes_n_Stage_II","genes_per_Stage_II","genes_padj_Stage_II","symbol_Stage_III","genes_n_Stage_III","genes_per_Stage_III","genes_padj_Stage_III")]


# Enriched tumor genes
enriched_stage_II<-rownames(expr_stage_II)[rownames(expr_stage_II) %in% pathways[["HALLMARK_MYC_TARGETS_V1"]]]
enriched_stage_III<-rownames(expr_stage_III)[rownames(expr_stage_III) %in% pathways[["HALLMARK_ANDROGEN_RESPONSE"]]]
enriched_stage_III<-rownames(expr_stage_III)[rownames(expr_stage_III) %in% pathways[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]]]

# Select genes_rankData_tumor
genes_rankData_tumor[genes_rankData_tumor$entrezgene_id %in% enriched_stage_II,"hgnc_symbol"]
genes_rankData_tumor[genes_rankData_tumor$entrezgene_id %in% enriched_stage_III,"hgnc_symbol"]
genes_rankData_tumor[genes_rankData_tumor$entrezgene_id %in% enriched_stage_III,"hgnc_symbol"]


##############################################################################################################

tumor_genes<-unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene))
tumor_samples<-c(sample_stage_I,sample_stage_II,sample_stage_III)

# Expression of tumor genes
expr_stage_tumor<-na.omit(df_reads_count_all_projects[[normalization_scheme]][tumor_genes,tumor_samples])

# Take for each ensembl_gene_id the entrezgene_accession, entrezgene_id, hgnc_symbol
genes_rankData_tumor     <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_tumor),mart=mart)

# Next, Take the annotation for each genes
genes_rankData_tumor   <-genes_rankData_tumor[genes_rankData_tumor$ensembl_gene_id %in% rownames(expr_stage_tumor),]

# Take the first occurance of each ensembl_gene_id 
genes_rankData_tumor   <- genes_rankData_tumor[match(unique(genes_rankData_tumor$ensembl_gene_id), genes_rankData_tumor$ensembl_gene_id),]

# Set the rownames as the ensembl_gene_id
rownames(genes_rankData_tumor)<-genes_rankData_tumor$ensembl_gene_id

# Set the rownames entrezgene_id
rownames(genes_rankData_tumor)   <- genes_rankData_tumor[rownames(genes_rankData_tumor),"entrezgene_id"]

# Set the rownames entrezgene_id
rownames(expr_stage_tumor)   <- genes_rankData_tumor[rownames(genes_rankData_tumor),"entrezgene_id"]

# Compute the geseca
geseca_Stage_tumor   <- data.frame(geseca(pathways, expr_stage_tumor))

# Filter padj
geseca_Stage_tumor  <-geseca_Stage_tumor[geseca_Stage_tumor$padj<0.05,]

# Enriched tumor genes
enriched_tumor_genes<-rownames(expr_stage_tumor)[rownames(expr_stage_tumor) %in% pathways[[geseca_Stage_tumor$pathway]]]

# Select genes_rankData_tumor
genes_rankData_tumor[genes_rankData_tumor$entrezgene_id %in% enriched_tumor_genes,"hgnc_symbol"]

###############################################################################################################################3
# Save TSV file with genes from Stage3
write_tsv(results_hallmark, paste(output_dir,"/hallmarks_genes.tsv",sep=""))

write_tsv(geseca_Stage_I, paste(output_dir,"/hallmarks_genes_stage_I.tsv",sep=""))
write_tsv(geseca_Stage_II, paste(output_dir,"/hallmarks_genes_stage_II.tsv",sep=""))
write_tsv(geseca_Stage_III, paste(output_dir,"/hallmarks_genes_stage_III.tsv",sep=""))


# check the hallmarks against the paper

# Perform gene set enrichment analysis
# Install msigdb packages.
# fgsea
# Genes from stage I     (numGene, n)    : number of genes of this stage associated to Hallmark X.
# Genes from stage I     (%)             : percentage of genes of this stage associated to Hallmark X.
# Genes from stage II    (numGene, n)    : number of genes of this stage associated to Hallmark X.
# Genes from stage II    (%)             : percentage of genes of this stage associated to Hallmark X.
# Genes from stage III   (numGene, n)    : number of genes of this stage associated to Hallmark X.
# Genes from stage III   (%)             : percentage of genes of this stage associated to Hallmark X.
# Stage III minus I      (Δ)             : number of genes enriched in stage III minus number of genes enriched in stage I
