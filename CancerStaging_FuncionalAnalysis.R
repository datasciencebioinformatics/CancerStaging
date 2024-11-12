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
write_tsv(all_anotation_results, paste(output_dir,"/hallmarks_genes.tsv",sep=""))

