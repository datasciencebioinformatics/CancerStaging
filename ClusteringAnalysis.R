# Add a collumns with the sample id
transporse_normalized_table$sample_id<-rownames(transporse_normalized_table)

# Merge transporse_normalized_table with metadatamerged_data_patient_info
transporse_normalized_table<-unique(merge(transporse_normalized_table,merged_data_patient_info,by="sample_id")[,c(genes,"tissue_type","sample_id")])

# set rownames
rownames(transporse_normalized_table)<-transporse_normalized_table$sample_id

# Split into trainning and testing
trainning <- sample_n(transporse_normalized_table,dim(transporse_normalized_table)[1]/2)
testing   <- transporse_normalized_table[which(!rownames(transporse_normalized_table) %in% rownames(trainning)),]

# Basic Parameter Tuning
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

# model 3 logistic regression
mod3 <- train(tissue_type ~ ., method = 'glm', data = trainning, trControl = fitControl, preProcess = c('center', 'scale'))

