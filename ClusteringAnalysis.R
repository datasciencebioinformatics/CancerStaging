# Add a collumns with the sample id
transporse_normalized_table$sample_id<-rownames(transporse_normalized_table)

# Merge transporse_normalized_table with metadata
transporse_normalized_table<-merge(transporse_normalized_table,merged_data_patient_info,by="sample_id")

transporse_normalized_table<-transporse_normalized_table[patient,c(genes,"tissue_type")]

# model 3 logistic regression
mod3 <- train(tissue_type ~ ., method = 'glm', data = train, trControl = tCtrl, preProcess = c('center', 'scale'))

# predict with mod 3
pred3 <- predict(mod3, newdata = test)

# confusionMatrix with mod 3
cfm3 <- confusionMatrix(pred3, test$booking_status)

# first roc curve
roc1 <- roc(test$booking_status, pred1Probability)

# calculate auc
auc(roc1)

# plot roc curves
plot(roc1, col = 'black', lty = 2, main = "ROC", legacy.axes = TRUE, percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", print.auc = TRUE)

# Plot_raw_vibration_data.png               
png(filename=paste(output_dir,"Plot_Roc_curve.png",sep=""), width = 20, height = 30, res=600, units = "cm")  
plot(roc1, col = 'black', lty = 2, main = "ROC", legacy.axes = TRUE, percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", print.auc = TRUE)
dev.off()
