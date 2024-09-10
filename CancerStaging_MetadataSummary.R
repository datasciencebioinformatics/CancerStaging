# Take the name of all varibles (all collumns)
variable_names<-colnames(merged_data_patient_info)

# Take the name of each sample (uniquely)
sample_id_names<-unique(merged_data_patient_info$sample_id)

# df_completeness
df_completeness<-data.frame(data.frame(Variable=c(),completeness=c()))

# for each variable
for (variable in variable_names)
{
  # Take all values of this variable
  variable_values<-unique(merged_data_patient_info[merged_data_patient_info$sample_id %in% sample_id_names,c(variable,"sample_id")])[,variable]

  # Number of valkues
  length_variables<-length(unique(merged_data_patient_info[merged_data_patient_info$sample_id %in% sample_id_names,c(variable,"sample_id")])[,variable])

  # Trim the values  
  variable_values<-str_trim(variable_values)

  # Remove empty space
  variable_values<-variable_values[variable_values != ""]

  # Remove NA space
  variable_values<-variable_values[variable_values != "NA"]  

  # Remove NA space
  variable_values<-variable_values[variable_values != "-"]    
  
  # Remove NA space
  variable_values<-variable_values[!is.na(variable_values)]

  # df_completeness
  df_completeness <- rbind(df_completeness,data.frame(Variable=variable,completeness=length(variable_values)/length_variables*100))
}
# df_completeness_only_copmplete
df_completeness_only_copmplete<-df_completeness[df_completeness$completeness>80,]

# df_completeness_only_copmplete
df_completeness_only_copmplete<-df_completeness_only_copmplete[which(!df_completeness_only_copmplete$Variable %in% c("Sample.Type","Sample.ID","Case.ID","Project.ID","Data.Type","Data.Category","File.Name","File.ID","case_submitter_id","project_id","sample_submitter_id","case_id","project_id.x","case_submitter_id.x","sample_type","sample_type_id","sample_id")),]


# Simple scatter plot
# FindClusters_resolution
png(filename=paste(output_dir,"df_completeness_only_copmplete.png",sep=""), width = 28, height = 14, res=600, units = "cm")
  ggplot(data=df_completeness_only_copmplete, aes(x=Variable, y=completeness)) + geom_bar(stat = "identity")
dev.off()


