####################################################################################################################
# A script to normalize reads count to RPKM                                                                        #
####################################################################################################################
# Gene length with TxDb.Hsapiens.UCSC.hg19.knownGene and getGeneLengthAndGCContent.                                #
# Plot correlation of length among the two.                                                                        #
####################################################################################################################
# Compute gene length with xDb.Hsapiens.UCSC.hg19.knownGene                                                        #
                                                                         # 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene                                                                          #
tx_by_gene <- transcriptsBy(txdb, by="gene")                                                                       #
gene_lens <- max(width(tx_by_gene))                                                                                #
                                                                                                                   #
# Compute gene getGeneLengthAndGCContent                                                                           #
getGeneLengthAndGCContent("ENSG00000012048", "hsa")                                                                #
####################################################################################################################


####################################################################################################################
# RPKM normalization                                                                                               #
# The normalization is done for each 1000 genes duo to limitation in the biomart connection                        #
# First, gene length and gc content for all genes in the reads count table                                         #
# Take the gene names, without variant identification                                                              #
# vector to store all gene ids                                                                                     #
df_gene_ids<-data.frame(gene_id=c(),gene_id_cp=c())                                                                #
for (gene_id in rownames(unstranded_data))                                                                         #
{                                                                                                                  #
    # Store gene ids                                                                                               #
    print(gene_id)
    gene_ids<-strsplit(gene_id,".",fixed=T)[[1]][[1]]                                                              #                                                  
                                                                                                                   #
    # Contatenate gene lists                                                                                       #
    df_gene_ids<-rbind(df_gene_ids,data.frame(gene_id=gene_ids,gene_id_cp=gene_id))                                #
}                                                                                                                  #
# Split gene_ids vector in parts                                                                                   #
gene_ids_vector<-split(df_gene_ids$gene_id,ceiling(seq_along(df_gene_ids$gene_id) / 1000))                         #
####################################################################################################################
# Data.frame to store geneLengthAndGCContent                                                                       #
df_geneLengthAndGCContent<-data.frame(length=c(),gc=c())                                                           #
                                                                                                                   #
# For each part of the vectors                                                                                     #
for (index in names(gene_ids_vector) )                                                                             #
{                                                                                                                  #
    # Concatenate files                                                                                            ########
    df_geneLengthAndGCContent<-rbind(df_geneLengthAndGCContent,getGeneLengthAndGCContent(gene_ids_vector[[index]], "hsa"))#
}                                                                                                                         #
rownames(df_geneLengthAndGCContent)[!grepl(".", rownames(df_geneLengthAndGCContent), fixed=TRUE)]                         #
###########################################################################################################################





# edgeR rpkm
y <- DGEList(counts=counts,genes=data.frame(Length=GeneLength))
y <- calcNormFactors(y)
RPKM <- rpkm(y)

# deseq2 fpkm
dds <- DESeqDataSetFromMatrix(countData = cts, colData = NULL, design= ~ 1)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
# the FPKM values
fpkm(dds)
