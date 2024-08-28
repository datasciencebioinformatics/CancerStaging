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
df_reads_count_per_project_fpkm$IDS<-substring(rownames(df_reads_count_per_project_fpkm),1,last=15)
df_reads_count_per_project_tpm$IDS<-substring(rownames(df_reads_count_per_project_tpm),1,last=15)

# Keep only the first occurance
df_reads_count_all_projects_raw <- reads_count_all_projects[match(unique(df_reads_count_all_projects_raw$IDS), df_reads_count_all_projects_raw$IDS),]
df_reads_count_per_project_fpkm <- reads_count_all_projects[match(unique(df_reads_count_per_project_fpkm$IDS), df_reads_count_per_project_fpkm$IDS),]
df_reads_count_per_project_tpm <- reads_count_all_projects[match(unique(df_reads_count_per_project_tpm$IDS), df_reads_count_per_project_tpm$IDS),]

# Rename cols
rownames(df_reads_count_all_projects_raw)<-df_reads_count_all_projects_raw$IDS
rownames(df_reads_count_per_project_fpkm)<-df_reads_count_per_project_fpkm$IDS
rownames(df_reads_count_per_project_tpm)<-df_reads_count_per_project_tpm$IDS

# Read count all project
df_reads_count_all_projects_raw<-df_reads_count_all_projects_raw[,!colnames(df_reads_count_all_projects_raw) %in% c("IDS")]
df_reads_count_per_project_fpkm<-df_reads_count_per_project_fpkm[,!colnames(df_reads_count_per_project_fpkm) %in% c("IDS")]
df_reads_count_per_project_tpm<-df_reads_count_per_project_tpm[,!colnames(df_reads_count_per_project_tpm) %in% c("IDS")]

# Use only genes for which gene length is available
df_reads_count_all_projects_raw<-df_reads_count_all_projects_raw[which(rownames(df_reads_count_all_projects_raw) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]
df_reads_count_per_project_fpkm<-df_reads_count_per_project_fpkm[which(rownames(df_reads_count_per_project_fpkm) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]
df_reads_count_per_project_tpm<-df_reads_count_per_project_tpm[which(rownames(df_reads_count_per_project_tpm) %in% geneLength_ENTREZID_ENSEMBL$ENSEMBL),]

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
unstranded_rpkm<-edgeR::rpkm(reads_count_all_projects[geneLength_ENTREZID_ENSEMBL$ENSEMBL,], gene.length = geneLength_ENTREZID_ENSEMBL$geneLength) #
##############################################################################################################
unstranded_dgelist <- DGEList(counts=reads_count_all_projects[geneLength_ENTREZID_ENSEMBL$ENSEMBL,],genes=data.frame(Length=geneLength_ENTREZID_ENSEMBL$geneLength))
unstranded_dgelist <- calcNormFactors(unstranded_dgelist, method = c("TMM")
unstranded_dgelist_rpkm <- edgeR::rpkm(unstranded_dgelist)
####################################################################################################################
# NOISeq
unstranded_NOISeq_rpkm = NOISeq::rpkm(reads_count_all_projects[geneLength_ENTREZID_ENSEMBL$ENSEMBL,], long = geneLength_ENTREZID_ENSEMBL$geneLength, lc = 1, k = 0)
####################################################################################################################
# TMM normalization with no length correction
unstranded_NOISeq_TMM = NOISeq::tmm(reads_count_all_projects[geneLength_ENTREZID_ENSEMBL$ENSEMBL,], long = geneLength_ENTREZID_ENSEMBL$geneLength, lc = 0, k = 0)
####################################################################################################################
## RPKM normalization
TP53_edgeR_rpkm_dgelist_TMM_TPM<-data.frame(RPKM_edgeR=unstranded_rpkm["ENSG00000141510",],RPKM_dgelist=unstranded_dgelist_rpkm["ENSG00000141510",], RPKM_NOISeq=unstranded_NOISeq_rpkm["ENSG00000141510",], TMM_from_NoiseSeq=unstranded_NOISeq_TMM["ENSG00000141510",])
####################################################################################################################           
# Calclation matrix
cor_TP53_RPKM_TMM_TPM <- cor(TP53_edgeR_rpkm_dgelist_TMM_TPM)
round(cor_TP53_RPKM_TMM_TPM, 2)
####################################################################################################################
# Save normalized data                                                                                             #
write_tsv(data.frame(unstranded_rpkm),         paste(output_dir,"unstranded_edgeR_rpkm.tsv",sep=""))			   #
write_tsv(data.frame(unstranded_dgelist_rpkm), paste(output_dir,"unstranded_dgelist_rpkm.tsv",sep=""))			   #
write_tsv(data.frame(unstranded_NOISeq_rpkm),   paste(output_dir,"unstranded_NOISeq_rpkm.tsv",sep=""))			   #
write_tsv(data.frame(unstranded_NOISeq_TMM),   paste(output_dir,"unstranded_NOISeq_TMM.tsv",sep=""))               #
####################################################################################################################

