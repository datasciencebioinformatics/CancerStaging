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
expr_stage_I<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_I_gene,sample_stage_I])
expr_stage_II<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_II_gene,sample_stage_II])
expr_stage_III<-na.omit(df_reads_count_all_projects[[normalization_scheme]][selected_genes_Stage_III_gene,sample_stage_III])

expr_stage_I<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_I,sample_stage_I])
expr_stage_II<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_II,sample_stage_II])
expr_stage_III<-na.omit(df_reads_count_all_projects[[normalization_scheme]][unique_stage_III,sample_stage_III])

genes_rankData_stage_I     <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_I),mart=mart)
genes_rankData_stage_II    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_II),mart=mart)
genes_rankData_stage_III   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id","hgnc_symbol"),values=rownames(expr_stage_III),mart=mart)

genes_rankData_stage_I   <-genes_rankData_stage_I[genes_rankData_stage_I$ensembl_gene_id %in% rownames(expr_stage_I),]
genes_rankData_stage_II  <-genes_rankData_stage_II[genes_rankData_stage_II$ensembl_gene_id %in% rownames(expr_stage_II),]
genes_rankData_stage_III <-genes_rankData_stage_III[genes_rankData_stage_III$ensembl_gene_id %in% rownames(expr_stage_III),]

genes_rankData_stage_I   <- genes_rankData_stage_I[match(unique(genes_rankData_stage_I$ensembl_gene_id), genes_rankData_stage_I$ensembl_gene_id),]
genes_rankData_stage_II  <- genes_rankData_stage_II[match(unique(genes_rankData_stage_II$ensembl_gene_id), genes_rankData_stage_II$ensembl_gene_id),]
genes_rankData_stage_III <- genes_rankData_stage_III[match(unique(genes_rankData_stage_III$ensembl_gene_id), genes_rankData_stage_III$ensembl_gene_id),]

rownames(genes_rankData_stage_I)<-genes_rankData_stage_I$ensembl_gene_id
rownames(genes_rankData_stage_II)<-genes_rankData_stage_II$ensembl_gene_id
rownames(genes_rankData_stage_III)<-genes_rankData_stage_III$ensembl_gene_id

rownames(expr_stage_I)   <- genes_rankData_stage_I[rownames(expr_stage_I),"entrezgene_id"]
rownames(expr_stage_II)  <- genes_rankData_stage_II[rownames(expr_stage_II),"entrezgene_id"]
rownames(expr_stage_III) <- genes_rankData_stage_III[rownames(expr_stage_III),"entrezgene_id"]

symgbol_genes_rankData_stage_I   <- genes_rankData_stage_I
symgbol_genes_rankData_stage_II  <- genes_rankData_stage_II
symgbol_genes_rankData_stage_III <- genes_rankData_stage_III

rownames(expr_stage_I_symbol)   <- genes_rankData_stage_I[rownames(expr_stage_I),"hgnc_symbol"]
rownames(expr_stage_II_symbol)  <- genes_rankData_stage_II[rownames(expr_stage_II),"hgnc_symbol"]
rownames(expr_stage_III_symbol) <- genes_rankData_stage_III[rownames(expr_stage_III),"hgnc_symbol"]


geseca_Stage_I   <- data.frame(geseca(pathways, expr_stage_I))
geseca_Stage_II  <- data.frame(geseca(pathways, expr_stage_II))
geseca_Stage_III <- data.frame(geseca(pathways, expr_stage_III))

# Filter padj
geseca_Stage_I[geseca_Stage_I$padj<0.05,]
geseca_Stage_II[geseca_Stage_II$padj<0.05,]
geseca_Stage_III[geseca_Stage_III$padj<0.05,]

# For each hallmark
for (hallmarks in names(pathways))
{
  # Take the genes
  pathways[[hallmarks]]

  symbol_Stage_I<-paste(genes_rankData_stage_I[rownames(expr_stage_I) %in% pathways[[hallmarks]],"hgnc_symbol"],collapse=" , ")
  symbol_Stage_II<-paste(genes_rankData_stage_II[rownames(expr_stage_II) %in% pathways[[hallmarks]],"hgnc_symbol"],collapse=" , ")
  symbol_Stage_III<-genes_rankData_stage_III[rownames(expr_stage_III) %in% pathways[[hallmarks]],"hgnc_symbol"]

  # Take number of genes from this ptahway on stage I
  genes_Stage_I<-paste(rownames(expr_stage_I)[rownames(expr_stage_I) %in% pathways[[hallmarks]]],collapse=" , ")
  genes_Stage_II<-paste(rownames(expr_stage_II)[rownames(expr_stage_II) %in% pathways[[hallmarks]]],collapse=" , ")
  genes_Stage_III<-paste(rownames(expr_stage_III)[rownames(expr_stage_III) %in% pathways[[hallmarks]]],collapse=" , ")

  # Take number of genes from this ptahway on stage I
  genes_n_Stage_I<-sum(rownames(expr_stage_I) %in% pathways[[hallmarks]])
  genes_n_Stage_II<-sum(rownames(expr_stage_II) %in% pathways[[hallmarks]])
  genes_n_Stage_III<-sum(rownames(expr_stage_III) %in% pathways[[hallmarks]])

  print(data.frame(hallmark=hallmarks,symbol_Stage_I=symbol_Stage_I,symbol_Stage_II=symbol_Stage_II,symbol_Stage_III=hallmarks))
  print(data.frame(hallmark=hallmarks,genes_Stage_I=genes_Stage_I,genes_Stage_II=genes_Stage_II,genes_Stage_III=genes_Stage_III))
  print(data.frame(hallmark=hallmarks,genes_n_Stage_I=genes_n_Stage_I,genes_n_Stage_II=genes_n_Stage_II,genes_n_Stage_III=genes_n_Stage_III))
}



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
