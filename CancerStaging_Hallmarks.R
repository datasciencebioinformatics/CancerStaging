library("msigdb")
library("fgsea")

# Install msigdb packages.
# use the custom accessor to select a specific version of MSigDB
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')

# retrieeve the hallmarks gene sets
# 50 hallmarks
hallmarks_gene_set<-subsetCollection(msigdb.hs, 'h')

# log2FC of genes for stages, I, II, III
rankData_stage_I   <- 
rankData_stage_II  <- 
rankData_stage_III <- 

# check the hallmarks against the paper
fgseaRes_stage_I   <- fgsea(hallmarks_gene_set, rankData_stage_I, minSize = 15, maxSize = 500)
fgseaRes_stage_II  <- fgsea(hallmarks_gene_set, rankData_stage_II, minSize = 15, maxSize = 500)
fgseaRes_stage_III <- fgsea(hallmarks_gene_set, rankData_stage_III, minSize = 15, maxSize = 500)


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
