makeFunctionalPrediction = function(
  path_to_otu_table,
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "Tax4Fun2_prediction",
  database_mode = 'Ref100NR',
  normalize_by_copy_number = T,
  min_identity_to_reference = 97,
  include_user_data = F,
  path_to_user_data = '',
  name_of_user_data = '',
  use_uproc = T,
  normalize_pathways = F,
  write_counts = F)
{
  # Read old log file to control the database mode
  path_to_log_file = file.path(path_to_temp_folder, 'logfile1.txt')
  log_file = read.delim(path_to_log_file, header = F)
  if(database_mode != as.character(log_file$V1[2])) stop('Logfile indicates other database mode')
  if(as.character(log_file$V1[4]) == "User data will be included" & include_user_data == FALSE) warning("Your log file indicates that you used user data in the first step\nSet include_user_data to TRUE to incorporate!")
  if(as.character(log_file$V1[4]) != "User data will be included" & include_user_data == TRUE)
  {
    warning("Your log file indicates that you did not used user data in the first step\nSet include_user_data will be set to FALSE!")
    include_user_data = FALSE
  }
      
  # Set path to reference files
  if(database_mode == 'Ref99NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref99NR')
  } else if (database_mode == 'Ref100NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref100NR')
  } else {
    stop('Database mode unknown! valid choces are Ref99NR and Ref100NR')
  }

  if(min_identity_to_reference < 1) min_identity_to_reference = min_identity_to_reference * 100
  if(min_identity_to_reference < 90) warning("Minimum identity of less than 90% will likly results in inaccurate predictions!")
  message(paste0("Using minimum idenity cutoff of ", min_identity_to_reference, "% to nearest neighbor"))
  # Write to new log file
  path_to_log_file = file.path(path_to_temp_folder, 'logfile2.txt')
  write(x = "Tax4fun2 beta", file = path_to_log_file, append = F)
  write(x = date(), file = path_to_log_file, append = T)

  # Reading and reducing the blast file
  ref_blast_result = read.delim(file.path(path_to_temp_folder, 'ref_blast.txt'), h = F)
  ref_blast_result_reduced = ref_blast_result[which(ref_blast_result$V3 >= min_identity_to_reference), 1:2]

  if(include_user_data) {
    ref_blast_result_reduced_v1 = as.character(ref_blast_result_reduced$V1)
    ref_blast_result_reduced_v2 = as.character(ref_blast_result_reduced$V2)
    user_blast_result = read.delim(paste(path_to_temp_folder, '/user_blast.txt', sep = ''), h = F)
    user_blast_result_reduced = user_blast_result[which(user_blast_result$V3 >= min_identity_to_reference), 1:2]
    for(i in 1:nrow(user_blast_result_reduced))
    {
      user_id = as.character(user_blast_result_reduced$V1)[i]
      if(user_id %in% ref_blast_result_reduced_v1)
      {
        j = which(ref_blast_result_reduced$V1 == user_id)
        ref_blast_result_reduced_v2[j] = as.character(user_blast_result_reduced$V2)[i]
      } else {
        ref_blast_result_reduced_v1 = c(ref_blast_result_reduced_v1, user_id)
        ref_blast_result_reduced_v2 = c(ref_blast_result_reduced_v2, as.character(user_blast_result_reduced$V2)[i])
      }
    }
    ref_blast_result_reduced = data.frame(V1 = ref_blast_result_reduced_v1, V2 = ref_blast_result_reduced_v2)
  }
  ref_blast_result_reduced
  
  # Reading and filtering the otu table
  otu_table = read.delim(path_to_otu_table)
  otu_table_reduced = merge(x = ref_blast_result_reduced, y = otu_table, by.x = 'V1', by.y = names(otu_table)[1])[,-1]
  otu_table_reduced_aggregated = aggregate(x = otu_table_reduced[,-1], by = list(otu_table_reduced[,1]), sum)

  # Write unknown fraction to log file
  if((ncol(otu_table) - 1) == 1)
  {
    unknown_fraction1 = as.data.frame(round(1 - sum(ifelse(otu_table_reduced[,-1]>0,1,0)) / sum(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - sum(otu_table_reduced_aggregated[,-1]) / sum(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  } else {
    unknown_fraction1 = as.data.frame(round(1 - colSums(ifelse(otu_table_reduced[,-1]>0,1,0)) / colSums(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - colSums(otu_table_reduced_aggregated[,-1]) / colSums(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  }

  # normalize or not normalize is the question
  n = 1
  if(use_uproc) n = 3
  if(normalize_by_copy_number) n = n + 1
  
  # Generate reference profile
  message('Generating reference profile')
  reference_profile = NULL
  for(reference_id in otu_table_reduced_aggregated$Group.1)
  {
    #print(reference_id)
    if(grepl(pattern = 'user', x = reference_id))
    {
      path_to_profile = file.path(path_to_user_data, name_of_user_data)
      reference_id = gsub("_[0-9]*", "", reference_id)
      reference_file_path = file.path(path_to_profile, paste0(reference_id, '.tbl'))
    } else{
      path_to_profile = path_to_ref_profiles
      reference_file_path = file.path(path_to_profile, paste0(reference_id, '.tbl.gz'))
    }
    reference_file = read.delim(file = reference_file_path)
    reference_profile = rbind(reference_profile, as.numeric(reference_file[,n]))
  }
  dim(reference_profile)

  ko_list = read.delim(file.path(path_to_reference_data, 'KEGG/ko.txt'))
  
  # Calculate functional profiles sample-wise
  message('Generating functional profile for:')
  functional_prediction = NULL
  functional_prediction_counts = NULL
  for(sample in 2:ncol(otu_table_reduced_aggregated))
  {
    message(names(otu_table_reduced_aggregated[sample]))
    functional_prediction_sample = reference_profile * as.numeric(otu_table_reduced_aggregated[,sample])
    functional_prediction_sample = colMeans(functional_prediction_sample)
    # Functional counts
    if(write_counts) functional_prediction_counts = cbind(functional_prediction_counts, functional_prediction_sample)
    # Functional percents
    functional_prediction_sample = functional_prediction_sample / sum(functional_prediction_sample)
    if(is.na(sum(functional_prediction_sample))) functional_prediction_sample[1:nrow(ko_list)] = 0
    functional_prediction = cbind(functional_prediction, functional_prediction_sample)
  }
  colnames(functional_prediction) = names(otu_table)[2:ncol(otu_table_reduced_aggregated)]
  if(write_counts) colnames(functional_prediction_counts) = names(otu_table)[2:ncol(otu_table_reduced_aggregated)]
  functional_prediction_final = data.frame(KO = ko_list$ko, functional_prediction, description = ko_list$description)
  if(write_counts) functional_prediction_counts_final = data.frame(KO = ko_list$ko, functional_prediction_counts, description = ko_list$description)
  if(ncol(functional_prediction) >= 2) keep = which(rowSums(functional_prediction) > 0)
  if(ncol(functional_prediction) == 1) keep = which(functional_prediction > 0)
  if (length(keep) == 0) stop("No functional prediction possible!\nEither no nearest neighbor found or your table is empty!")
  functional_prediction_final = functional_prediction_final[keep,]
  if(write_counts) functional_prediction_counts_final = functional_prediction_counts_final[keep,]
  write.table(x = functional_prediction_final, file = file.path(path_to_temp_folder, 'functional_prediction.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  if(write_counts) write.table(x = functional_prediction_counts_final, file = file.path(path_to_temp_folder, 'functional_prediction_counts.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # Converting the KO profile to a profile of KEGG pathways
  message('Converting functions to pathways')
  if(write_counts) functional_prediction = functional_prediction_counts
  ko2ptw = read.delim(file.path(path_to_reference_data, 'KEGG/ko2ptw.txt'))
  if(normalize_pathways) functional_prediction = functional_prediction / ko_list$ptw_count
  pathway_prediction = aggregate(x = functional_prediction[ko2ptw$nrow,], by = list(ko2ptw$ptw), sum)
  if(write_counts) pathway_prediction_counts = pathway_prediction
  if(ncol(pathway_prediction) >= 3)
  {
    col_sums = colSums(pathway_prediction[,-1])
    col_sums[col_sums == 0] = 1
    pathway_prediction[,-1] = t(t(pathway_prediction[,-1]) / col_sums)
    keep = which(rowSums(pathway_prediction[,-1]) > 0)
  } else {
    pathway_prediction[,-1] = t(t(pathway_prediction[,-1]) / sum(pathway_prediction[,-1]))
    keep = which(pathway_prediction[,2] > 0)
  }
  if(sum(pathway_prediction[,-1]) == 0) stop("Conversion to pathway failed!")
  names(pathway_prediction) = names(otu_table)
  names(pathway_prediction)[1] = 'pathway'
  if(write_counts) names(pathway_prediction_counts) = names(pathway_prediction)

  ptw_desc = read.delim(paste(path_to_reference_data, '/KEGG/ptw.txt', sep = ''))
  pathway_prediction_final = merge(pathway_prediction, ptw_desc)[keep,]
  if(write_counts) pathway_prediction_counts_final = merge(pathway_prediction_counts, ptw_desc)[keep,]
  write.table(x = pathway_prediction_final, file = file.path(path_to_temp_folder, 'pathway_prediction.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  if(write_counts) write.table(x = pathway_prediction_counts_final, file = file.path(path_to_temp_folder, 'pathway_prediction_counts.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
}
