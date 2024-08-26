m####################################################################################################################
# A script to normalize reads count to RPKM                                                                        #
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
    print(gene_id)                                                                                                 #
    gene_ids<-strsplit(gene_id,".",fixed=T)[[1]][[1]]                                                              #                                                  
                                                                                                                   #
    # Contatenate gene lists                                                                                       #
    df_gene_ids<-rbind(df_gene_ids,data.frame(gene_id=gene_ids,gene_id_cp=gene_id))                                #
}                                                                                                                  #
###########################################################################################################################
# Here, use the goseq package to retrieve the gene length
getlength_vector<-getlength(df_gene_ids$gene_id,'hg19','ensGene')

# getlength_df
getlength_df<-data.frame(ENSEMBL=df_gene_ids$gene_id,getlength=getlength_vector)
####################################################################################################################                                                                       # 
# Gene length with TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx_by_gene <- transcriptsBy(txdb, by="gene")
gene_lens <- max(width(tx_by_gene))

# Take the gene length
txdb_geneLength<-data.frame(ENTREZID=names(gene_lens), geneLength=as.vector(gene_lens))

# ids_stage_I - all ENSEMBL anotated using bitr
gene_emsembl      <-bitr(txdb_geneLength$ENTREZID, fromType = "ENTREZID", toType = c("ENSEMBL","SYMBOL"), OrgDb="org.Hs.eg.db")

# Merge tables
txdb_geneLength<-merge(txdb_geneLength,gene_emsembl,by="ENTREZID")
####################################################################################################################
merge_geneLengh_Counts<-merge(getlength_df,txdb_geneLength,by="ENSEMBL")
####################################################################################################################
# Remove NA lines
merge_geneLengh_Counts<-na.omit(merge_geneLengh_Counts)

cor(merge_geneLengh_Counts$getlength,merge_geneLengh_Counts$geneLength)







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
