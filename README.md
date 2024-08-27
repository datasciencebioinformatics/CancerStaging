A cohort with the following TCGA projects was constructed from the https://portal.gdc.cancer.gov: 

Breast/Breast Invasive Carcinoma/TCGA-BRCA
Lung/Lung Adenocarcinoma/TCGA-LUAD
Lung/Lung Squamous Cell Carcinoma/TCGA-LUSC
Colon/Colon Adenocarcinoma/TCGA-COAD/
Rectum/Rectum Adenocarcinoma TCGA-READ/
Prostate/Prostate Adenocarcinoma TCGA-PRAD/
Skin/Skin Cutaneous Melanoma/TCGA-SKCM/
Stomach/Stomach Adenocarcinoma/TCGA-STAD/
Liver/Liver Hepatocellular Carcinoma/TCGA-LIHC/

These cancer types were selected for having greatest incidence by 2020, acording to WHO: https://www.who.int/news-room/fact-sheets/detail/cancer. Sub--types of lung cancer LUSC and LUAD were selected because our interest and comparing sub-types and different tissue.

#### 0- Pre-configuration mount rclone to access google cloud server
rclone config                                       

rclone mount googledrive: /home/felipe/googledrive/ 

#### 1- Download tables from gdc_manifest.2024-07-18.txt
/home/felipe/Documents/Cancer_staging/github/CancerStaging_DownloadFromManifest.sh

#### 2- Create expression tables from file
/home/felipe/Documents/Cancer_staging/github/CancerStaging_CreateTableFromFiles.sh

#### 3- CancerStaging Setup All Paramters
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")

#### 4- Load all R packages
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")

#### 5- Create metadata CancerStaging_CreateMetadataFromGDCFiles
##### The following files were downloaded from GDC portal, clinical.tsv (information for each patient), sample.tsv (information for each sample) and exposure.tsv (data about patient life-style). These three tables were merged by the "case_id" value. Additionally, gdc_sample_sheet.2024-07-18.tsv were also downloaded (infomation about the experimental files). Then, the gdc_sample_sheet.2024-07-18.tsv was merge with the formed merged table, relating the "Sample ID" to "sample_submitter_id", respectivelly. 
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_CreateMetadataFromGDCFiles.R")

#### 6- Load expression data CancerStaging_LoadExpressionData
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadExpressionData.R")

#### 7- Load expression data CancerStaging_LoadExpressionData. Raw read counts are normalized with DESeq2 package for CPM (counts per million), TPM (transcripts per kilobase million) and RPKM/FPKM (reads/fragments per kilobase of exon per million reads/fragments mapped).
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_ExpressionDataNormalization.R")

#### 8- Create Paired Samples Tumor Normal
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_CreatePairedSamplesTumorNormal.R")

#### 9- Load interactome data
##### The IntAct interactome was obtained from the intact-micluster.txt file (version updated December 2017) accessed on January 11, 2018, with 152280 interactions among 15651 gene symbols. After converting gene symbols to ENSEMBL identifiers with EnsemblToUniprotKBconversionList.txt, 148169 interactions (97.3%) and 14492 genes (92.6%) were kept. To calculate the connectivity per gene, we counted the number of times each gene appeared in the interactome. 
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadInteractomeData.R")

#### 10- Create Expression maps ExpressionMaps
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_ExpressionMaps.R")


