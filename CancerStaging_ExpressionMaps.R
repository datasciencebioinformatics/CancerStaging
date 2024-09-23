# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_scheme<-"tpm"

# Construnction of 3d coordinates.
# X=T2           : Thymine composition in second codon position (T2)
# Y=connectivity : interactome_data                 # Carels checked    .
# Z=expression   : reads_count_all_projects         # Use the fpkm
# Reading the contents of TSV file using read_tsv() method
Interactomes_GC3_T2_file <-"/home/felipe/Documents/github/CancerStaging/Interactomes_GC3_T2.csv"

# Read data
Interactomes_GC3_T2_data <-read.table(file = Interactomes_GC3_T2_file, sep = '\t', header = TRUE,fill=TRUE) 

# Rename second collumn
colnames(Interactomes_GC3_T2_data)[2]<-"SYMBOL"

# Merge Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL
Interactomes_GC3_T2_merged<-merge(geneLength_ENTREZID_ENSEMBL,Interactomes_GC3_T2_data,by="SYMBOL")

# Accross cancer types.
# The hypohtoses is to put the clusters ordered by entropy in the 3d maps                   .
# The hypotheseis is to order the clusters by entropy in the 3d map                         .
# Amont the 2d coordinates this can have low signal, in a three map can have a stroger sinal.
expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

# I need take the intersection of ids from the interactome table and expression table.
# The expression table has unique ids but it seems the Interactomes_GC3_T2_merged is duplicated (check).
# Uniprot entries have multiple ENSEMBL. Criteria is to take UNIRPOTKB with highest conections.
# To Do : take UNIRPOTKB with highest conections.

# I want this table today by 16:00.
# Input for the table. For each gene:
# 1) Expression per patient. 2) T2, Connection, etc.


# Set ensembl ids
ENSEMBL_ids<-unique(intersect(rownames(expression_table_normalized),Interactomes_GC3_T2_merged$ENSEMBL))
                                 
# Take the used genes
used_genes<-unique(rownames(expression_table_normalized[which(rownames(expression_table_normalized) %in% Interactomes_GC3_T2_merged$ENSEMBL),]))

Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$ENSEMBL %in% used_genes,]

dim(expression_table_normalized[used_genes,])





[Interactomes_GC3_T2_merged$ENSEMBL,]



Interactomes_GC3_T2_merged<-merge(expression_table_normalized,Interactomes_GC3_T2_data,by="SYMBOL")



