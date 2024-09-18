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
# Path to CDS files
cds_file <- "/home/felipe/Downloads/cds.txt"                                                  

# Load data
cds_data <-read.table(file = cds_file, sep = ' ', header = TRUE,fill=TRUE)    
####################################################################################################################
# Then I will merge the dataset to obtain the ENSEMBL ids
id_conversion_ENTREZID_ENSEMBL <-bitr(store_gene_length$ENTREZID, fromType = "ENTREZID", toType = c("SYMBOL","ENSEMBL","UNIPROT"), OrgDb="org.Hs.eg.db")
####################################################################################################################
geneLength_ENTREZID_ENSEMBL<-merge(store_gene_length,id_conversion_ENTREZID_ENSEMBL,by="ENTREZID")
####################################################################################################################
geneLength_ENTREZID_ENSEMBL<-merge(cds_data,geneLength_ENTREZID_ENSEMBL,by="UNIPROT")
colnames(geneLength_ENTREZID_ENSEMBL)[2]<-"geneLength"
geneLength_ENTREZID_ENSEMBL<-geneLength_ENTREZID_ENSEMBL[,c("ENTREZID","geneLength","SYMBOL","ENSEMBL")]
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
#############################################################################################################################
df_reads_count_all_projects_rpkm<-data.frame(edgeR::rpkm(df_reads_count_all_projects_raw[rownames(geneLength_ENTREZID_ENSEMBL),], gene.length = geneLength_ENTREZID_ENSEMBL$geneLength)) #
#############################################################################################################################
# Calculate RPK
# RPK = RAW RNA COUNTS / GENE LENGTH
RPK <- df_reads_count_all_projects_raw[rownames(geneLength_ENTREZID_ENSEMBL),] / geneLength_ENTREZID_ENSEMBL$geneLength

# Replace infinite by NA
# Infinite or can happen if gene length of read counts are zero
RPK[sapply(RPK, is.infinite)] <- NA

# Infinite and NA values are removed to compute the sum of RPK values over the genessssss
df_reads_count_all_projects_tpm_calculated <- t( t(RPK) * 1e6 / colSums(RPK, na.rm=TRUE) )
##########################################################################################################################
colnames(df_reads_count_all_projects_tmm)<-colnames(df_reads_count_all_projects_raw)
colnames(df_reads_count_all_projects_rpkm)<-colnames(df_reads_count_all_projects_raw)
colnames(df_reads_count_all_projects_tpm_calculated)<-colnames(df_reads_count_all_projects_tpm_calculated)                                                      

colnames(df_reads_count_all_projects_tmm)<-colnames(df_reads_count_all_projects_raw)
colnames(df_reads_count_all_projects_rpkm)<-colnames(df_reads_count_all_projects_raw)
colnames(df_reads_count_all_projects_tpm_calculated)<-colnames(df_reads_count_all_projects_tpm_calculated)                                                      
colnames(df_reads_count_all_projects_tpm)<-colnames(df_reads_count_all_projects_raw)                                                      
colnames(df_reads_count_all_projects_rpkm)<-colnames(df_reads_count_all_projects_raw)                                                      
##########################################################################################################################
rownames(df_reads_count_all_projects_tmm)<-rownames(df_reads_count_all_projects_raw)
rownames(df_reads_count_all_projects_rpkm)<-rownames(df_reads_count_all_projects_raw)
rownames(df_reads_count_all_projects_tpm_calculated)<-rownames(df_reads_count_all_projects_tpm_calculated)                                                      
rownames(df_reads_count_all_projects_tpm)<-rownames(df_reads_count_all_projects_raw)                                                      
rownames(df_reads_count_all_projects_rpkm)<-rownames(df_reads_count_all_projects_raw)              
##########################################################################################################################
#write.table(df_reads_count_all_projects_tmm,              paste(output_dir,"df_reads_count_all_projects_tmm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#write.table(data.frame(df_reads_count_all_projects_rpkm), paste(output_dir,"df_reads_count_all_projects_rpkm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
##########################################################################################################################
tp53_raw<-t(data.frame(df_reads_count_all_projects_raw["ENSG00000141510",]))
tp53_tmm<-t(data.frame(df_reads_count_all_projects_tmm["ENSG00000141510",]))
tp53_fpkm<-t(data.frame(df_reads_count_all_projects_fpkm["ENSG00000141510",]))
tp53_tpm<-data.frame(df_reads_count_all_projects_tpm["ENSG00000141510",])
tp53_rpkm<-t(data.frame(df_reads_count_all_projects_rpkm["ENSG00000141510",]))
tp53_tpm_calc<-data.frame(df_reads_count_all_projects_tpm_calculated["ENSG00000141510",])

colnames(tp53_tpm)<-colnames(tp53_raw)
colnames(tp53_tpm_calc)<-colnames(tp53_raw)

df_normalization<-data.frame(raw=tp53_raw[,"ENSG00000141510"],tp53_tmm[,"ENSG00000141510"],tp53_fpkm[,"ENSG00000141510"],tp53_tpm[,"ENSG00000141510"],tp53_rpkm[,"ENSG00000141510"],tp53_tpm_calc[,"ENSG00000141510"])
###################################################################################
colnames(df_normalization)<-c("raw","tmm","fpkm","tpm","rpkm","tpm_calc")
####################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"df_normalization.png",sep=""), width = 24, height = 24, res=600, units = "cm")
  ggpairs(df_normalization,sgnf=3)
dev.off()
####################################################################################
df_reads_count_all_projects<-list(raw=df_reads_count_all_projects_raw,tmm=df_reads_count_all_projects_tmm,fpkm=df_reads_count_all_projects_fpkm,tpm=df_reads_count_all_projects_tpm,rpkm=df_reads_count_all_projects_rpkm, tpm_calc=df_reads_count_all_projects_tpm_calculated)
print("\nCancerStaging_ExpressionDataNormalization")
