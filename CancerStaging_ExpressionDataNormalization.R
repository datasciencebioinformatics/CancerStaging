####################################################################################################################
# A script to normalize reads count to RPKM                                                                         #
# TxDb.Hsapiens.UCSC.hg19.knownGene is used to obtain the length of CCDS sequence data. Data for 26649 transcripts are present in the database.
# There are 60660 entries in the gene set,and 60616 unique trancripts.
####################################################################################################################
# Gene length with TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Take the transcript Lengths by the with.cds_len
Hsapiens_length <- transcriptLengths(txdb, with.cds_len=TRUE,with.utr5_len=TRUE,with.utr3_len=TRUE)

# Take the gene length
txdb_geneLength<-data.frame(ENTREZID=Hsapiens_length$gene_id, geneLength=Hsapiens_length$cds_len)

# ids_stage_I - all ENSEMBL anotated using bitr
gene_emsembl      <-bitr(txdb_geneLength$ENTREZID, fromType = "ENTREZID", toType = c("ENSEMBL","SYMBOL"), OrgDb="org.Hs.eg.db")

# Merge tables
txdb_geneLength<-merge(txdb_geneLength,gene_emsembl,by="ENTREZID")
####################################################################################################################
# store_gene_length
store_gene_length<-data.frame(ENTREZID=c(),geneLength=c())

# for each ENTREZID in txdb_geneLength$ENTREZID
for (ENTREZID in unique(txdb_geneLength$ENTREZID))
{
    # Take the ENTREZID
    ENTREZID_ID  <-unique(txdb_geneLength[which(txdb_geneLength$ENTREZID %in% ENTREZID),"ENTREZID"])
   
    # Take the maximum gene length
    geneLength<-max(txdb_geneLength[which(txdb_geneLength$ENTREZID %in% ENTREZID_ID),"geneLength"])

    # Store gene length
    store_gene_length<-rbind(store_gene_length,data.frame(ENTREZID=ENTREZID_ID,geneLength=geneLength))
}
####################################################################################################################
# Then I will merge the dataset to obtain the ENSEMBL ids
id_conversion_ENTREZID_ENSEMBL <-bitr(store_gene_length$ENTREZID, fromType = "ENTREZID", toType = c("SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")
####################################################################################################################
geneLength_ENTREZID_ENSEMBL<-merge(store_gene_length,id_conversion_ENTREZID_ENSEMBL,by="ENTREZID")
####################################################################################################################
# Substrintg of dataset
df_reads_count_all_projects_raw$IDS<-substring(rownames(df_reads_count_all_projects_raw),1,last=15)
df_reads_count_all_projects_fpkm$IDS<-substring(rownames(df_reads_count_all_projects_fpkm),1,last=15)
df_reads_count_all_projects_tpm$IDS<-substring(rownames(df_reads_count_all_projects_tpm),1,last=15)

# Keep only the first occurance
df_reads_count_all_projects_raw <- df_reads_count_all_projects_raw[match(unique(df_reads_count_all_projects_raw$IDS), df_reads_count_all_projects_raw$IDS),]
df_reads_count_all_projects_fpkm <- df_reads_count_all_projects_fpkm[match(unique(df_reads_count_all_projects_fpkm$IDS), df_reads_count_all_projects_fpkm$IDS),]
df_reads_count_all_projects_tpm <- df_reads_count_all_projects_tpm[match(unique(df_reads_count_all_projects_tpm$IDS), df_reads_count_all_projects_tpm$IDS),]

# Rename cols
rownames(df_reads_count_all_projects_raw)<-df_reads_count_all_projects_raw$IDS
rownames(df_reads_count_all_projects_fpkm)<-df_reads_count_all_projects_fpkm$IDS
rownames(df_reads_count_all_projects_tpm)<-df_reads_count_all_projects_tpm$IDS

# Read count all project
df_reads_count_all_projects_raw<-df_reads_count_all_projects_raw[,!colnames(df_reads_count_all_projects_raw) %in% c("IDS")]
df_reads_count_all_projects_fpkm<-df_reads_count_all_projects_fpkm[,!colnames(df_reads_count_all_projects_fpkm) %in% c("IDS")]
df_reads_count_all_projects_tpm<-df_reads_count_all_projects_tpm[,!colnames(df_reads_count_all_projects_tpm) %in% c("IDS")]

# Use only genes for which gene length is available
df_reads_count_all_projects_raw<-df_reads_count_all_projects_raw[which(rownames(df_reads_count_all_projects_raw) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]
df_reads_count_all_projects_fpkm<-df_reads_count_all_projects_fpkm[which(rownames(df_reads_count_all_projects_fpkm) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]
df_reads_count_all_projects_tpm<-df_reads_count_all_projects_tpm[which(rownames(df_reads_count_all_projects_tpm) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]

####################################################################################################################
#Sort gene length data.frame
# Keep only the first occurance
geneLength_ENTREZID_ENSEMBL <- geneLength_ENTREZID_ENSEMBL[match(unique(geneLength_ENTREZID_ENSEMBL$ENSEMBL), geneLength_ENTREZID_ENSEMBL$ENSEMBL),]

# First set row names
rownames(geneLength_ENTREZID_ENSEMBL)<-geneLength_ENTREZID_ENSEMBL$ENSEMBL

# Now, check the row names
rownames(geneLength_ENTREZID_ENSEMBL)<-geneLength_ENTREZID_ENSEMBL$ENSEMBL

# Remove NA entries
geneLength_ENTREZID_ENSEMBL<-geneLength_ENTREZID_ENSEMBL[which(!is.na(geneLength_ENTREZID_ENSEMBL$ENSEMBL)),]

# Set rownames
rownames(geneLength_ENTREZID_ENSEMBL)<-geneLength_ENTREZID_ENSEMBL$ENSEMBL

# Sort table according to count table
geneLength_ENTREZID_ENSEMBL<-geneLength_ENTREZID_ENSEMBL[rownames(df_reads_count_all_projects_raw),]
####################################################################################################################
unstranded_dgelist              <- DGEList(counts=df_reads_count_all_projects_raw[geneLength_ENTREZID_ENSEMBL$ENSEMBL,],genes=data.frame(Length=geneLength_ENTREZID_ENSEMBL$geneLength))
unstranded_dgelist              <- calcNormFactors(unstranded_dgelist, method = c("TMM"))
df_reads_count_all_projects_tmm <- data.frame(cpm(unstranded_dgelist))
###########################################################################################################################
write_tsv(df_reads_count_all_projects_tmm, "/home/felipe/Documents/Cancer_staging/df_reads_count_all_projects_tmm.tsv") #
##########################################################################################################################
df_tp53_expressuion<-data.frame(raw=df_reads_count_all_projects_raw["ENSG00000141510",],tmm=df_reads_count_all_projects_tmm["ENSG00000141510",], fpkm=df_reads_count_all_projects_fpkm["ENSG00000141510",], tmm=df_reads_count_all_projects_tmm["ENSG00000141510",])

tp53_raw<-data.frame(df_reads_count_all_projects_raw["ENSG00000141510",])
tp53_tmm<-data.frame(df_reads_count_all_projects_tmm["ENSG00000141510",])
tp53_fpkm<-data.frame(df_reads_count_all_projects_tmm["ENSG00000141510",])
tp53_tpm<-data.frame(df_reads_count_all_projects_tmm["ENSG00000141510",])

# Data.frames with gene name and read counts
df_raw <-data.frame(as.vector(df_reads_count_all_projects_raw))
df_raw$gene<-rownames(df_reads_count_all_projects_raw)

# Data.frames with gene name and read counts
df_fpkm <-data.frame(as.vector(df_reads_count_all_projects_fpkm))
df_fpkm$gene<-rownames(df_reads_count_all_projects_fpkm)

# Data.frames with gene name and read counts
df_tmm <-data.frame(as.vector(df_reads_count_all_projects_tmm))
df_tmm$gene<-rownames(df_reads_count_all_projects_tmm)

# Data.frames with gene name and read counts
df_tpm <-data.frame(as.vector(df_reads_count_all_projects_tpm))
df_tpm$gene<-rownames(df_reads_count_all_projects_tpm)
allall_nr

merge(merge(df_fpkm,df_tpm,by="gene"),merge(df_tpm,df_tmm,by="gene"),by="gene")
                   
merge(df_raw,merge(df_fpkm,by="gene"),)

3
data.frame(gene=data.frame(gene=rownames(df_reads_count_all_projects_tpm),tpm,=as.vector(df_reads_count_all_projects_tpm)))
data.frame(gene=data.frame(gene=rownames(df_reads_count_all_projects_tmm),tmm=as.vector(df_reads_count_all_projects_tmm)))

cbind(, df_reads_count_all_projects_tmm["ENSG00000141510",], df_reads_count_all_projects_fpkm["ENSG00000141510",], df_reads_count_all_projects_tmm["ENSG00000141510",])

pairs(~ raw + tmm + fpkm + tmm, data = as.matrix(df_tp53_expressuion))





