####################################################################################################################
# cds_file <- "/home/felipe/Downloads/cds.txt"                                                  
# Load data
# cds_data <-read.table(file = cds_file, sep = ' ', header = TRUE,fill=TRUE)    
####################################################################################################################
# Then I will merge the dataset to obtain the ENSEMBL ids
# id_conversion_ENTREZID_ENSEMBL <-bitr(store_gene_length$ENTREZID, fromType = "ENTREZID", toType = c("SYMBOL","ENSEMBL","UNIPROT"), OrgDb="org.Hs.eg.db")
####################################################################################################################
# geneLength_ENTREZID_ENSEMBL<-merge(store_gene_length,id_conversion_ENTREZID_ENSEMBL,by="ENTREZID")
####################################################################################################################
# geneLength_ENTREZID_ENSEMBL<-merge(cds_data,geneLength_ENTREZID_ENSEMBL,by="UNIPROT")
# geneLength_ENTREZID_ENSEMBL<-geneLength_ENTREZID_ENSEMBL[,c("ENTREZID","geneLength","SYMBOL","ENSEMBL")]
####################################################################################################################
