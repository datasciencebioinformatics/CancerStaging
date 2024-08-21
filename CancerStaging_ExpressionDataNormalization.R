####################################################################################################################
# A script to normalize reads count to RPKM                                                                        #
####################################################################################################################
# Gene length with TxDb.Hsapiens.UCSC.hg19.knownGene and getGeneLengthAndGCContent.                                #
# Plot correlation of length among the two.                                                                        #
####################################################################################################################
# Compute gene length with xDb.Hsapiens.UCSC.hg19.knownGene                                                        #
library(TxDb.Hsapiens.UCSC.hg19.knownGene)                                                                         # 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene                                                                          #
tx_by_gene <- transcriptsBy(txdb, by="gene")                                                                       #
gene_lens <- max(width(tx_by_gene))                                                                                #
                                                                                                                   #
# Compute gene getGeneLengthAndGCContent                                                                           #
getGeneLengthAndGCContent("ENSG00000012048", "hsa")                                                                #
####################################################################################################################
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
