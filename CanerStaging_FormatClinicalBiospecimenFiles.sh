# Organize bioscpecimen and clinical txt data
cp /home/felipe/Downloads/biospecimen.project-tcga-brca.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/
cp /home/felipe/Downloads/clinical.project-tcga-brca.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/

cp /home/felipe/Downloads/clinical.project-tcga-luad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/
cp /home/felipe/Downloads/biospecimen.project-tcga-luad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/

cp /home/felipe/Downloads/clinical.project-tcga-luad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/
cp /home/felipe/Downloads/biospecimen.project-tcga-luad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/

cp /home/felipe/Downloads/clinical.project-tcga-coad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/colon/Colon_Adenocarcinoma_TCGA_COAD/
cp /home/felipe/Downloads/biospecimen.project-tcga-coad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/colon/Colon_Adenocarcinoma_TCGA_COAD/

cp /home/felipe/Downloads/clinical.project-tcga-read.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/
cp /home/felipe/Downloads/biospecimen.project-tcga-read.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/

cp /home/felipe/Downloads/clinical.project-tcga-prad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/
cp /home/felipe/Downloads/biospecimen.project-tcga-prad.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/

cp /home/felipe/Downloads/clinical.project-tcga-skcm.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/
cp /home/felipe/Downloads/biospecimen.project-tcga-skcm.2024-07-16.tar.gz /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/

cp /home/felipe/Downloads/clinical.project-tcga-stad.2024-07-16.tar.gz  /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/	
cp /home/felipe/Downloads/biospecimen.project-tcga-stad.2024-07-16.tar.gz  /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/

cp /home/felipe/Downloads/clinical.project-tcga-lihc.2024-07-16.tar.gz  /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/
cp /home/felipe/Downloads/biospecimen.project-tcga-lihc.2024-07-16.tar.gz  /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/		

# Untar the zip files
tar xzf /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/biospecimen.project-tcga-brca.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/
tar xzf /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/clinical.project-tcga-brca.2024-07-16.tar.gz -C  /home/felipe/googledrive/Cancer_staging/breast/Breast_Invasive_Carcinoma_TCGA_BRCA/

tar xzf /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/clinical.project-tcga-luad.2024-07-16.tar.gz -C  /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/
tar xzf /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/biospecimen.project-tcga-luad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/lung/Lung_Adenocarcinoma_TCGA_LUAD/

tar xzf /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/clinical.project-tcga-luad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/
tar xzf /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/biospecimen.project-tcga-luad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/lung/Lung_Squamous_Cell_Carcinoma_TCGA_LUSC/

tar xzf /home/felipe/googledrive/Cancer_staging/colon//Colon_Adenocarcinoma_TCGA_COAD/clinical.project-tcga-coad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/colon/Colon_Adenocarcinoma_TCGA_COAD/
tar xzf /home/felipe/googledrive/Cancer_staging/colon/Colon_Adenocarcinoma_TCGA_COAD/biospecimen.project-tcga-coad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/colon/Colon_Adenocarcinoma_TCGA_COAD/

tar xzf /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/clinical.project-tcga-read.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/
tar xzf /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/biospecimen.project-tcga-read.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/rectum/Rectum_Adenocarcinoma_TCGA_READ/

tar xzf /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/clinical.project-tcga-prad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/
tar xzf /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/biospecimen.project-tcga-prad.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/prostate/Prostate_Adenocarcinoma_TCGA_PRAD/

tar xzf /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/clinical.project-tcga-skcm.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/
tar xzf /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/biospecimen.project-tcga-skcm.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/skin/Skin_Cutaneous_Melanoma_TCGA_SKCM/

tar xzf /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/clinical.project-tcga-stad.2024-07-16.tar.gz -C  /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/	
tar xzf /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/biospecimen.project-tcga-stad.2024-07-16.tar.gz -C  /home/felipe/googledrive/Cancer_staging/stomach/Stomach_Adenocarcinoma_TCGA_STAD/

tar xzf /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/clinical.project-tcga-lihc.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/
tar xzf /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/biospecimen.project-tcga-lihc.2024-07-16.tar.gz -C /home/felipe/googledrive/Cancer_staging/liver/Liver_Hepatocellular_Carcinoma_TCGA_LIHC/		
