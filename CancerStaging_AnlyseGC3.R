x<-c(20,50,70,90)

stage_I_genes<-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage I","GC3"]
stage_II_genes<-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage II","GC3"]
stage_III_genes<-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage III","GC3"]

ENSEMBL_stage_I  <-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage I","ENSEMBL"]
ENSEMBL_stage_II <-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage II","ENSEMBL"]
ENSEMBL_stage_III<-Interactomes_GC3_T2_values[Interactomes_GC3_T2_values$Stage=="Stage III","ENSEMBL"]

stage_I_classification<-data.frame(ENSEMBL=ENSEMBL_stage_I,class=cut(stage_I_genes, breaks=quantile(stage_I_genes, probs = c(0,0.2,0.5,0.7,0.9,1)), labels=c("0-20","20-50","50-70","70-90","90-100")))


stage_II_classification<-data.frame(ENSEMBL=ENSEMBL_stage_II,class=cut(stage_II_genes, breaks=quantile(stage_II_genes, probs = c(0,0.2,0.5,0.7,0.9,1)), labels=c("0-20","20-50","50-70","70-90","90-100")))


stage_III_classification<-data.frame(ENSEMBL=ENSEMBL_stage_III,class=cut(stage_III_genes, breaks=quantile(stage_III_genes, probs = c(0,0.2,0.5,0.7,0.9,1)), labels=c("0-20","20-50","50-70","70-90","90-100")))

stage_I_classification<-merge(Interactomes_GC3_T2_merged_Stage_I,stage_I_classification,by="ENSEMBL")
stage_II_classification<-merge(Interactomes_GC3_T2_merged_Stage_II,stage_II_classification,by="ENSEMBL")
stage_III_classification<-merge(Interactomes_GC3_T2_merged_Stage_III,stage_III_classification,by="ENSEMBL")


