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
reads_count_all_projects$IDS<-substring(rownames(reads_count_all_projects),1,last=15)

# Keep only the first occurance
reads_count_all_projects <- reads_count_all_projects[match(unique(reads_count_all_projects$IDS), reads_count_all_projects$IDS),]

# Rename cols
rownames(reads_count_all_projects)<-reads_count_all_projects$IDS

# Read count all project
reads_count_all_projects<-reads_count_all_projects[,!colnames(reads_count_all_projects) %in% c("IDS")]

# Use only genes for which gene length is available
reads_count_all_projects<-reads_count_all_projects[which(rownames(reads_count_all_projects) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]

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
geneLength_ENTREZID_ENSEMBL<-geneLength_ENTREZID_ENSEMBL[rownames(reads_count_all_projects),]

####################################################################################################################
# TO DO tomorrow 24th-august-2024
# RPKM normalization
# Check meticulously the DGEList normalization
# Important to check if the gene length use is correct.
# RPKM
unstranded_rpkm<-rpkm(reads_count_all_projects[geneLength_ENTREZID_ENSEMBL$ENSEMBL,], gene.length = geneLength_ENTREZID_ENSEMBL$geneLength) #
##############################################################################################################
# TPM normalization
reads_count_TPM     <- apply(reads_count_RPKM, 2, function(x) x / sum(as.numeric(x)) * 10^6) %>% as.data.frame()

# TPM
VeroTPM <- convertCounts(reads_count_all_projects, unit       = "tpm", geneLength = geneLength_ENTREZID_ENSEMBL,  log        = FALSE, prior.count=0)
####################################################################################################################
# TMM normalization
reads_count_DGEList <- calcNormFactors(reads_count_DGEList, method = "TMM")
reads_count_TMM <- cpm(reads_count_DGEList)
####################################################################################################################
# Normalizaton matrix
TP53_RPKM_TMM_TPM<-data.frame(RPKM=reads_count_RPKM["ENSG00000141510",],TPM=reads_count_TPM["ENSG00000141510",],TMM=reads_count_TMM["ENSG00000141510",])

# Calclation matrix
cor_TP53_RPKM_TMM_TPM <- cor(TP53_RPKM_TMM_TPM)
round(cor_TP53_RPKM_TMM_TPM, 2)
                         

