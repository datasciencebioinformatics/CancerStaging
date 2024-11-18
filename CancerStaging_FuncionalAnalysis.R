############################################################################################################################################
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# ENSEMBL ids were converted to ENTREZ ids. enrichGO on org.Hs.eg.db was used (pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, minGSSize = 3) to anotate 23, 62 and 169 genes from stages I, II and III, respectivelly. Then cnetplot was used to show asociations of genes to top 10 categories.
#######################################################################################################################################
genes_ids_stage_I<-unique_stage_I
genes_ids_stage_II<-unique_stage_II
genes_ids_stage_III<-unique_stage_III
########################################################################################################################################
# ids_stage_I - all ENSEMBL anotated using bitr
ids_stage_I      <-bitr(genes_ids_stage_I, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_II     <-bitr(genes_ids_stage_II, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_III    <-bitr(genes_ids_stage_III, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
########################################################################################################################################
# Set colnames,  gene_id, ENTREZID, SYMBOL
colnames(ids_stage_I)   <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_II)  <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_III) <-c("gene_id","ENTREZID","SYMBOL")
########################################################################################################################################
# Merge genes_unique and ids_stage in the same tables
#genes_Stage_I  <-merge(genes_unique_Stage_I,ids_stage_I,by="gene_id")
#genes_Stage_II <-merge(genes_unique_Stage_II,ids_stage_II,by="gene_id")
#genes_Stage_III<-merge(genes_unique_Stage_III,ids_stage_III,by="gene_id")

# Bind all three stages into one
genes_Stage_ALL<-unique(rbind(ids_stage_I,ids_stage_II,ids_stage_III))

# And set rownames(genes_Stage_ALL)
rownames(genes_Stage_ALL)<-genes_Stage_ALL$ENTREZID
########################################################################################################################################
# EnrichGO to obtain GO annotation, minGSSize = 3
go_ALL_Stage_I    <- enrichGO(gene = ids_stage_I$SYMBOL,  OrgDb  = org.Hs.eg.db,      ont = "ALL", pAdjustMethod = "BH",readable = TRUE,keyType = "SYMBOL", minGSSize = 3)
go_ALL_Stage_II   <- enrichGO(gene = ids_stage_II$SYMBOL, OrgDb  = org.Hs.eg.db,      ont = "ALL",  pAdjustMethod = "BH",readable = TRUE,keyType = "SYMBOL", minGSSize = 3)
go_ALL_Stage_III  <- enrichGO(gene = ids_stage_III$SYMBOL,OrgDb  = org.Hs.eg.db,      ont = "ALL",   pAdjustMethod = "BH",readable = TRUE,keyType = "SYMBOL", minGSSize = 3)
########################################################################################################################################
# enrichKEGG to obtain KEGG annotation, minGSSize = 3
kegg_ALL_Stage_I    <- enrichKEGG(gene = ids_stage_I$ENTREZ,  organism     = 'hsa', minGSSize = 3)
kegg_ALL_Stage_II   <- enrichKEGG(gene = ids_stage_II$ENTREZ,  organism     = 'hsa', minGSSize = 3)
kegg_ALL_Stage_III  <- enrichKEGG(gene = ids_stage_III$ENTREZ,  organism     = 'hsa', minGSSize = 3)
########################################################################################################################################
# enrichPathway to obtain enrichPathway annotation, minGSSize = 3
reactome_ALL_Stage_I    <- enrichPathway(gene=ids_stage_I$ENTREZ, readable=TRUE, minGSSize = 3)
reactome_ALL_Stage_II   <- enrichPathway(gene=ids_stage_II$ENTREZ, readable=TRUE, minGSSize = 3)
reactome_ALL_Stage_III  <- enrichPathway(gene=ids_stage_III$ENTREZ, readable=TRUE, minGSSize = 3)
########################################################################################################################################
# A table with the associate go term per gene.
# First, take the results per stage
# Go terms
go_results_Stage_I<-go_ALL_Stage_I@result
go_results_Stage_II<-go_ALL_Stage_II@result
go_results_Stage_III<-go_ALL_Stage_III@result

# Kegg terms
kegg_results_Stage_I<-kegg_ALL_Stage_I@result
kegg_results_Stage_II<-kegg_ALL_Stage_II@result
kegg_results_Stage_III<-kegg_ALL_Stage_III@result

# Reactome terms
reactome_results_Stage_I<-reactome_ALL_Stage_I@result
reactome_results_Stage_II<-reactome_ALL_Stage_II@result
reactome_results_Stage_III<-reactome_ALL_Stage_III@result

# Second, add information about stage
go_results_Stage_I$Stage<-"Stage I"
go_results_Stage_II$Stage<-"Stage II"
go_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
kegg_results_Stage_I$Stage<-"Stage I"
kegg_results_Stage_II$Stage<-"Stage II"
kegg_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
reactome_results_Stage_I$Stage<-"Stage I"
reactome_results_Stage_II$Stage<-"Stage II"
reactome_results_Stage_III$Stage<-"Stage III"

# Merge all stages
go_results_all_Stages         <-rbind(go_results_Stage_I,go_results_Stage_II,go_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]
kegg_results_all_Stages       <-rbind(kegg_results_Stage_I,kegg_results_Stage_II,kegg_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]
reactome_results_all_Stages   <-rbind(reactome_results_Stage_I,reactome_results_Stage_II,reactome_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]

# Merge all stages
go_results_all_Stages$Layer<-"GO"
kegg_results_all_Stages$Layer<-"KEGG"
reactome_results_all_Stages$Layer<-"Reactome"

# Merge all tables
all_anotation_results<-rbind(go_results_all_Stages,kegg_results_all_Stages,reactome_results_all_Stages)

# Filter by p.adjust
all_anotation_results<-all_anotation_results[which(all_anotation_results$p.adjust <0.05),]
######################################################################################################################
# Save TSV file with genes from Stage3
write_tsv(all_anotation_results, paste(output_dir,"/all_anotation_results.tsv",sep=""))

# Table wth
df_genes_terms<-data.frame(ID=c(),p.adjust=c(), Description=c(), geneID=c(), Count=c(),Stage=c(),Layer=c(),gene=c())

# For each term
for (term in unique(all_anotation_results$ID))
{
  # All the genes are take
  gene_IDs<-all_anotation_results[all_anotation_results$ID==term,"geneID"]

  ID            = all_anotation_results[all_anotation_results$ID==term,"ID"]
  p.adjust      = all_anotation_results[all_anotation_results$ID==term,"p.adjust"]
  Description   = all_anotation_results[all_anotation_results$ID==term,"Description"]
  geneID        = all_anotation_results[all_anotation_results$ID==term,"geneID"]
  Count         = all_anotation_results[all_anotation_results$ID==term,"Count"]
  Stage         = all_anotation_results[all_anotation_results$ID==term,"Stage"]
  Layer         = all_anotation_results[all_anotation_results$ID==term,"Layer"]  

  # For each gene, repeat the line
  for (gene in unlist(strsplit(gene_IDs,"/",fixed=T)) )
  {
    # Table wth
    df_genes_terms<-rbind(df_genes_terms,data.frame(ID=ID,p.adjust=p.adjust, Description=Description, geneID=geneID, Count=Count,Stage=Stage,Layer=Layer,gene=gene))
  }
}
######################################################################################################################
# Translate kegg back to symbols
# Convert ids
ids_stage_kegg      <-bitr(df_genes_terms[df_genes_terms$Layer=="KEGG","gene"], fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Match the ids
df_match_ids<-data.frame(ENTREZ=df_genes_terms[df_genes_terms$Layer=="KEGG","gene"],SYMOBL="")

# For each id
for (ENTREZ in df_match_ids$ENTREZ)
{
  print(ENTREZ)
  df_match_ids[which(df_match_ids$ENTREZ==ENTREZ),"SYMOBL"]<-ids_stage_kegg[ids_stage_kegg$ENTREZ==ENTREZ,"SYMBOL"]    
}
# I must match ids_stage_kegg with 
df_genes_terms[df_genes_terms$Layer=="KEGG","gene"]<-df_match_ids$SYMOBL
######################################################################################################################
# Take the 3 most abundat layers from each stage
# First, spearate per stage
df_genes_terms_stage_I   <-df_genes_terms[df_genes_terms$Stage=="Stage I",]
df_genes_terms_stage_II  <-df_genes_terms[df_genes_terms$Stage=="Stage II",]
df_genes_terms_stage_III <-df_genes_terms[df_genes_terms$Stage=="Stage III",]

# Second, separate by 
tabble_terms_GO_Stage_I       <-table(df_genes_terms_stage_I[df_genes_terms_stage_I$Layer=="GO","ID"])
tabble_terms_KEGG_Stage_I     <-table(df_genes_terms_stage_I[df_genes_terms_stage_I$Layer=="KEGG","ID"])
tabble_terms_Reactome_Stage_I <-table(df_genes_terms_stage_I[df_genes_terms_stage_I$Layer=="Reactome","ID"])

# Second, separate by 
tabble_terms_GO_Stage_II       <-table(df_genes_terms_stage_II[df_genes_terms_stage_II$Layer=="GO","ID"])
tabble_terms_KEGG_Stage_II     <-table(df_genes_terms_stage_II[df_genes_terms_stage_II$Layer=="KEGG","ID"])
tabble_terms_Reactome_Stage_II <-table(df_genes_terms_stage_II[df_genes_terms_stage_II$Layer=="Reactome","ID"])

# Second, separate by 
tabble_terms_GO_Stage_III       <-table(df_genes_terms_stage_III[df_genes_terms_stage_III$Layer=="GO","ID"])
tabble_terms_KEGG_Stage_III     <-table(df_genes_terms_stage_III[df_genes_terms_stage_III$Layer=="KEGG","ID"])
tabble_terms_Reactome_Stage_III <-table(df_genes_terms_stage_III[df_genes_terms_stage_III$Layer=="Reactome","ID"])

tabble_terms_GO_Stage_I<-tabble_terms_GO_Stage_I[order(-tabble_terms_GO_Stage_I)]
tabble_terms_Reactome_Stage_I<-tabble_terms_Reactome_Stage_I[order(-tabble_terms_Reactome_Stage_I)]
tabble_terms_KEGG_Stage_I<-tabble_terms_KEGG_Stage_I[order(-tabble_terms_KEGG_Stage_I)]

tabble_terms_GO_Stage_II<-tabble_terms_GO_Stage_II[order(-tabble_terms_GO_Stage_II)]
tabble_terms_Reactome_Stage_II<-tabble_terms_Reactome_Stage_II[order(-tabble_terms_Reactome_Stage_II)]
tabble_terms_KEGG_Stage_II<-tabble_terms_KEGG_Stage_II[order(-tabble_terms_KEGG_Stage_II)]

tabble_terms_GO_Stage_III<-tabble_terms_GO_Stage_III[order(-tabble_terms_GO_Stage_III)]
tabble_terms_Reactome_Stage_III<-tabble_terms_Reactome_Stage_III[order(-tabble_terms_Reactome_Stage_III)]
tabble_terms_KEGG_Stage_III<-tabble_terms_KEGG_Stage_III[order(-tabble_terms_KEGG_Stage_III)]

# Order the three most abundant terms by padd
tabble_terms_GO_Stage_I<-data.frame(tabble_terms_GO_Stage_I)
tabble_terms_GO_Stage_II<-data.frame(tabble_terms_GO_Stage_II)
tabble_terms_GO_Stage_III<-data.frame(tabble_terms_GO_Stage_III)

tabble_terms_KEGG_Stage_I<-data.frame(tabble_terms_KEGG_Stage_I)
tabble_terms_KEGG_Stage_II<-data.frame(tabble_terms_KEGG_Stage_II)
tabble_terms_KEGG_Stage_III<-data.frame(tabble_terms_KEGG_Stage_III)

tabble_terms_Reactome_Stage_I<-data.frame(tabble_terms_Reactome_Stage_I)
tabble_terms_Reactome_Stage_II<-data.frame(tabble_terms_Reactome_Stage_II)
tabble_terms_Reactome_Stage_III<-data.frame(tabble_terms_Reactome_Stage_III)

tabble_terms_GO_Stage_I$Padj<-all_anotation_results[tabble_terms_GO_Stage_I$Var1,"p.adjust"]
tabble_terms_GO_Stage_II$Padj<-all_anotation_results[tabble_terms_GO_Stage_II$Var1,"p.adjust"]
tabble_terms_GO_Stage_III$Padj<-all_anotation_results[tabble_terms_GO_Stage_III$Var1,"p.adjust"]

tabble_terms_KEGG_Stage_I$Padj<-all_anotation_results[tabble_terms_KEGG_Stage_I$Var1,"p.adjust"]
tabble_terms_KEGG_Stage_II$Padj<-all_anotation_results[tabble_terms_KEGG_Stage_II$Var1,"p.adjust"]
tabble_terms_KEGG_Stage_III$Padj<-all_anotation_results[tabble_terms_KEGG_Stage_III$Var1,"p.adjust"]

tabble_terms_Reactome_Stage_I$Padj<-all_anotation_results[tabble_terms_Reactome_Stage_I$Var1,"p.adjust"]
tabble_terms_Reactome_Stage_II$Padj<-all_anotation_results[tabble_terms_Reactome_Stage_II$Var1,"p.adjust"]
tabble_terms_Reactome_Stage_III$Padj<-all_anotation_results[tabble_terms_Reactome_Stage_III$Var1,"p.adjust"]


######################################################################################################################
# Merge all layers
tabble_terms_all_Stage_I<-rbind(rbind(tabble_terms_GO_Stage_I,tabble_terms_Reactome_Stage_I))
tabble_terms_all_Stage_II<-rbind(tabble_terms_GO_Stage_II,tabble_terms_KEGG_Stage_II)
tabble_terms_all_Stage_III<-rbind(tabble_terms_GO_Stage_III,tabble_terms_Reactome_Stage_III)

# Merge all layers
tabble_terms_all_Stage_I$Stage<-"Stage I"
tabble_terms_all_Stage_II$Stage<-"Stage II"
tabble_terms_all_Stage_III$Stage<-"Stage III"

# Merge all stages
table_terms_all_Stages<-rbind(tabble_terms_all_Stage_I,tabble_terms_all_Stage_II,tabble_terms_all_Stage_III)

# Add list of genes
table_terms_all_Stages$Description<-""
table_terms_all_Stages$Genes<-""

# For each term, take the genes
for (term in table_terms_all_Stages$Var1)
{
  # Store genes
  table_terms_all_Stages[table_terms_all_Stages$Var1==term,"Description"]<-unique(df_genes_terms[df_genes_terms$ID==term,"Description"])
  
  # Store genes
  table_terms_all_Stages[table_terms_all_Stages$Var1==term,"Genes"]<-unique(paste(df_genes_terms[df_genes_terms$ID==term,"gene"],collapse=", "))
}
# Save TSV file with genes from Stage3
write_tsv(table_terms_all_Stages, paste(output_dir,"/table_terms_all_terms_Stages.tsv",sep=""))
######################################################################################################################
######################################################################################################################
tabble_terms_GO_Stage_I<-na.omit(tabble_terms_GO_Stage_I[order(-tabble_terms_GO_Stage_I$Freq),][1:10,])
tabble_terms_GO_Stage_II<-na.omit(tabble_terms_GO_Stage_II[order(-tabble_terms_GO_Stage_II$Freq),][1:10,])
tabble_terms_GO_Stage_III<-na.omit(tabble_terms_GO_Stage_III[order(-tabble_terms_GO_Stage_III$Freq),][1:10,])

#tabble_terms_KEGG_Stage_I<-na.omit(tabble_terms_KEGG_Stage_I[order(-tabble_terms_KEGG_Stage_I$Freq),][1:10,])
tabble_terms_KEGG_Stage_II<-na.omit(tabble_terms_KEGG_Stage_II[order(-tabble_terms_KEGG_Stage_II$Freq),][1:10,])
#tabble_terms_KEGG_Stage_III<-na.omit(tabble_terms_KEGG_Stage_III[order(-tabble_terms_KEGG_Stage_III$Freq),][1:10,])

tabble_terms_Reactome_Stage_I<-na.omit(tabble_terms_Reactome_Stage_I[order(-tabble_terms_Reactome_Stage_I$Freq),][1:10,])
#tabble_terms_Reactome_Stage_II<-na.omit(tabble_terms_Reactome_Stage_II[order(-tabble_terms_Reactome_Stage_II$Freq),][1:10,])
#tabble_terms_Reactome_Stage_III<-na.omit(tabble_terms_Reactome_Stage_III[order(-tabble_terms_Reactome_Stage_III$Freq),][1:10,])
######################################################################################################################
# Merge all layers
tabble_terms_all_Stage_I<-rbind(rbind(tabble_terms_GO_Stage_I,tabble_terms_Reactome_Stage_I))
tabble_terms_all_Stage_II<-rbind(tabble_terms_GO_Stage_II,tabble_terms_KEGG_Stage_II)
tabble_terms_all_Stage_III<-rbind(tabble_terms_GO_Stage_III,tabble_terms_Reactome_Stage_III)

# Merge all layers
tabble_terms_all_Stage_I$Stage<-"Stage I"
tabble_terms_all_Stage_II$Stage<-"Stage II"
tabble_terms_all_Stage_III$Stage<-"Stage III"

# Merge all stages
table_terms_all_Stages<-rbind(tabble_terms_all_Stage_I,tabble_terms_all_Stage_II,tabble_terms_all_Stage_III)

# Add list of genes
table_terms_all_Stages$Description<-""
table_terms_all_Stages$Genes<-""

# For each term, take the genes
for (term in table_terms_all_Stages$Var1)
{
  # Store genes
  table_terms_all_Stages[table_terms_all_Stages$Var1==term,"Description"]<-unique(df_genes_terms[df_genes_terms$ID==term,"Description"])
  
  # Store genes
  table_terms_all_Stages[table_terms_all_Stages$Var1==term,"Genes"]<-unique(paste(df_genes_terms[df_genes_terms$ID==term,"gene"],collapse=", "))
}
# Save TSV file with genes from Stage3
write_tsv(table_terms_all_Stages, paste(output_dir,"/table_terms_all_Stages.tsv",sep=""))
######################################################################################################################
# To DO : 19-November2024
# Organize this table to see differences
# Colllums to represent each stage
######################################################################################################################
# Take all terms from STAGE I
stage_I_Terms<-df_genes_terms[df_genes_terms$Stage=="Stage I","ID"]
stage_II_Terms<-df_genes_terms[df_genes_terms$Stage=="Stage II","ID"]
stage_III_Terms<-df_genes_terms[df_genes_terms$Stage=="Stage III","ID"]

# Take the insersection from the three stages

# Take the most significants of the first stage
# Take the most significants of the second stage
# Take the most significants of the third stage
