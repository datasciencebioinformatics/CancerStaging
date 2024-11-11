library("msigdb")
library("fgsea")

# Load mart tables
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Install msigdb packages.
# use the custom accessor to select a specific version of MSigDB
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')

# retrieeve the hallmarks gene sets
# 50 hallmarks
hallmarks_gene_set<-subsetCollection(msigdb.hs, 'h')

# log2FC of genes for stages, I, II, III
df_rankData_stage_I   <- na.omit(df_FC[selected_genes_Stage_I_gene,c("Gene","FC")])
df_rankData_stage_II  <- na.omit(df_FC[selected_genes_Stage_II_gene,c("Gene","FC")])
df_rankData_stage_III <- na.omit(df_FC[selected_genes_Stage_III_gene,c("Gene","FC")]) 

vc_rankData_stage_I   <- df_rankData_stage_I$FC
vc_rankData_stage_II  <- df_rankData_stage_II$FC
vc_rankData_stage_III <- df_rankData_stage_III$FC

names(vc_rankData_stage_I)   <- df_rankData_stage_I$Gene
names(vc_rankData_stage_II)  <- df_rankData_stage_II$Gene
names(vc_rankData_stage_III) <- df_rankData_stage_III$Gene

genes_rankData_stage_I    <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id"),values=names(vc_rankData_stage_I),mart=mart)
genes_rankData_stage_II   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id"),values=names(vc_rankData_stage_II),mart=mart)
genes_rankData_stage_III  <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_accession","entrezgene_id"),values=names(vc_rankData_stage_III),mart=mart)

genes_rankData_stage_I   <-genes_rankData_stage_I[genes_rankData_stage_I$ensembl_gene_id %in% names(vc_rankData_stage_I),]
genes_rankData_stage_II  <-genes_rankData_stage_II[genes_rankData_stage_II$ensembl_gene_id %in% names(vc_rankData_stage_II),]
genes_rankData_stage_III <-genes_rankData_stage_III[genes_rankData_stage_III$ensembl_gene_id %in% names(vc_rankData_stage_III),]

rownames(genes_rankData_stage_I)<-genes_rankData_stage_I$ensembl_gene_id
rownames(genes_rankData_stage_II)<-genes_rankData_stage_II$ensembl_gene_id
rownames(genes_rankData_stage_III)<-genes_rankData_stage_III$ensembl_gene_id

genes_rankData_stage_I[names(vc_rankData_stage_I)]


# check the hallmarks against the paper
fgseaRes_stage_I   <- fgsea(hallmarks_gene_set, vc_rankData_stage_I, minSize = 15, maxSize = 500)
fgseaRes_stage_II  <- fgsea(hallmarks_gene_set, vc_rankData_stage_II, minSize = 15, maxSize = 500)
fgseaRes_stage_III <- fgsea(hallmarks_gene_set, vc_rankData_stage_III, minSize = 15, maxSize = 500)

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
