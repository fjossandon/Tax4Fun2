runRefBlast = function(path_to_otus, path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2_prediction", database_mode = 'Ref100NR', use_force = F, num_threads = 1, include_user_data = F, path_to_user_data = '', name_of_user_data = "")
{
  blast_bin = "blastn"
  res = system(command = paste(blast_bin, "-help"), intern = T)
  if(length(res) == 0)
  {
    blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn")
    if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn.exe")
    res = system(command = paste(blast_bin, "-help"), intern = T)
    if(length(res) == 0) stop("blastn not found! Consider to use the buildDependencies() command.")
    system(paste("chmod +x", blast_bin))
  }
  
  makeblastdb_bin = "makeblastdb"
  res = system(command = paste(blast_bin, "-help"), intern = T)
  if(length(res) == 0)
  {
    blast_bin = file.path(path_to_reference_data, "blast_bin/bin/makeblastdb")
    if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/bin/makeblastdb.exe")
    res = system(command = paste(blast_bin, "-help"), intern = T)
    if(length(res) == 0) stop("blastn not found! Consider to use the buildDependencies() command.")
    system(paste("chmod +x", makeblastdb_bin))
  }

  # Function to reduce blast results
  blastTableReducer = function(path_to_blast_file = '') {
    if(file.size(path_to_blast_file) == 0) stop('Blast file empty!')
    
    id1 = ""
    file_in = file(description = path_to_blast_file, open = "r")
    file_out = file(description = paste0(path_to_blast_file, ".tmp"), open = "w")
    while (TRUE) 
    {
      line = readLines(con = file_in, n = 1)
      if (length(line) == 0) break
      
      id2 = strsplit(x = line, split = "\t", fixed = T)[[1]][1]
      if(id1 != id2)
      {
        id1 = id2
        write(x = line, file = file_out, append = T)
      }
    }
    close(file_in)
    close(file_out)
    
    # Remove old file and rename tmp file
    file.remove(path_to_blast_file)
    file.rename(from = paste0(path_to_blast_file, ".tmp"), to = path_to_blast_file)
  }

  # Check for the presence of the otu file and the reference data
  if(!file.exists(path_to_otus)) stop('Otu file not found!')
  if(!dir.exists(path_to_reference_data)) stop('Reference data not found!')

  # Choose which refernence data set is used
  if(database_mode == 'Ref99NR') {
    path_to_ref_db = file.path(path_to_reference_data, 'Ref99NR/Ref99NR.fasta')
  } else if (database_mode == 'Ref100NR') {
    path_to_ref_db = file.path(path_to_reference_data, 'Ref100NR/Ref100NR.fasta')
  } else {
    stop('Database mode unknown! valid choces are Ref99NR, Ref100NR')
  }
  if(!file.exists(path_to_ref_db)) stop('Reference database not found!')
  
  # Check for the presence of the tmp folder. If present, stop or delete, if not, create
  if(dir.exists(path_to_temp_folder))
  {
    if(use_force) unlink(x = path_to_temp_folder, recursive = T)
    if(!use_force) stop('Temporay folder exist! Use \'use_force = T\' to overwrite')
  }
  dir.create(path_to_temp_folder)

  # Write to log file
  path_to_log_file = file.path(path_to_temp_folder, 'logfile1.txt')
  write(x = "RefBlast", file = path_to_log_file, append = F)
  write(x = database_mode, file = path_to_log_file, append = T)
  write(x = path_to_otus, file = path_to_log_file, append = T)
  if(include_user_data) write(x = 'User data will be included', file = path_to_log_file, append = T)
  write(x = date(), file = path_to_log_file, append = T)

  message('Copy and generate database')

  cmd = paste(makeblastdb_bin, '-dbtype nucl -in', path_to_ref_db)
  if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
  message('Reference blast started')
  cmd = paste(blast_bin, '-db', path_to_ref_db, '-query', path_to_otus, '-evalue 1e-20 -max_target_seqs 1000000 -outfmt 6 -out', file.path(path_to_temp_folder, 'ref_blast.txt'), '-num_threads', num_threads)
  if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
  blastTableReducer(file.path(path_to_temp_folder, 'ref_blast.txt'))
  message('Reference blast finished')
  message('Cleanup finished')

  if(include_user_data)
  {
    message('Copy and generate user database')
    path_to_user_db = file.path(path_to_user_data, name_of_user_data, 'USER.fasta')
    file.copy(from = path_to_user_db, to = file.path(path_to_temp_folder, 'user_blast.fasta'))
    cmd = paste(makeblastdb_bin, '-dbtype nucl -in', file.path(path_to_temp_folder, 'user_blast.fasta'))
    if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
    if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
    message('User blast started')
    cmd = paste(blast_bin, '-db', file.path(path_to_temp_folder, '/user_blast.fasta'), '-query', path_to_otus, '-evalue 1e-20 -max_target_seqs 1000000 -outfmt 6 -out', file.path(path_to_temp_folder, 'user_blast.txt'), '-num_threads', num_threads)
    if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, show.output.on.console = F)
    if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
    blastTableReducer(file.path(path_to_temp_folder, 'user_blast.txt'))
    message('User blast finished')
    unlink(file.path(path_to_temp_folder, 'user_blast.fasta*'))
    message('Cleanup finished')
  }
}
