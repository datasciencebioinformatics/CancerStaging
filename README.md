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
#### 1- Download tables from gdc_manifest.2024-07-18.txt
https://github.com/datasciencebioinformatics/CancerStaging/blob/main/CancerStaging_DownloadFromManifest.sh

#### 2- Create expression tables from file
https://github.com/datasciencebioinformatics/CancerStaging/blob/main/CancerStaging_CreateTableFromFiles.sh

#### 3- Load all R packages
source("https://github.com/datasciencebioinformatics/CancerStaging/blob/main/CancerStaging_LoadRPackages.R")

#### 4- Create metadata CancerStaging_CreateMetadataFromGDCFiles
source("https://github.com/datasciencebioinformatics/CancerStaging/blob/main/CancerStaging_CreateMetadataFromGDCFiles.R")

#### 5- Load expression data CancerStaging_LoadExpressionData
source("https://github.com/datasciencebioinformatics/CancerStaging/edit/main/CancerStaging_LoadExpressionData.R")
