setwd("E:/KY-INBRE PostDoc/Yvonne/Matt/RNAseq_Run1/Cell_Clustering")

Folder_list <- c("C2", "C3", "C4", "C7", "C12", "C13", "C14")

for (folder in Folder_list){
  files <- list.files(folder, pattern = "\\.txt$", full.name = TRUE )
  
  if (length(files) >0){
    merged_data <- do.call(rbind, lapply(files, read.table, header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 1))
    merged_data <- merged_data[!duplicated(merged_data[, 1]), ]
    
    # Split the content of the fourth column by "|" and keep the first part
    merged_data[, 4] <- sapply(strsplit(merged_data[, 4], "\\|"), `[`, 1)
    
    # Create a new data frame with the first and updated fourth column
    output_data <- merged_data[, c(1, 4)]
    
    # Save this data to a new file
    output_file <- file.path(folder, paste0(folder, "_merged_GeneSymbols.txt"))
    write.table(output_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
}

rawCounts <- read.table("E:/KY-INBRE PostDoc/Yvonne/Matt/RNAseq_Run1/COUNTS/htseq_rawCounts.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")

for (folder in Folder_list){
  if (folder == "C2"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR10", "MR11", "MR12", "MR13", "MR14", "MR15", "MR16", "MR17", "MR18")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR10", "MR11", "MR12")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR13", "MR14", "MR15")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR16", "MR17", "MR18")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C3"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR19", "MR20", "MR21", "MR22", "MR23", "MR24", "MR25", "MR27", "MR28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR19", "MR20", "MR21")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR22", "MR23", "MR24")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR25", "MR27", "MR28")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C4"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR28", "MR29", "MR30", "MR31", "MR32", "MR33")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR28", "MR29", "MR30")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR31", "MR32", "MR33")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C7"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR34", "MR35", "MR36", "MR37", "MR38", "MR39", "MR40", "MR41", "MR42")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR34", "MR35", "MR36")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR37", "MR38", "MR39")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR40", "MR41", "MR42")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C12"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR43", "MR44", "MR45", "MR46", "MR47", "MR48", "MR49", "MR50", "MR51")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR43", "MR44", "MR45")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR46", "MR47", "MR48")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR49", "MR50", "MR51")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C13"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR52", "MR53", "MR54", "MR55", "MR56", "MR57", "MR58", "MR59", "MR60")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR52", "MR53", "MR54")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR55", "MR56", "MR57")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR58", "MR59", "MR60")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  else if (folder == "C14"){
    readfile <- paste0("./", folder, "/", folder, "_merged_GeneSymbols.txt")
    df <- read.table(readfile, header = TRUE, sep = "\t")
    df_merged <- merge(df, rawCounts,  by.x = "V1", by.y = "ENSEMBL_ID", all.x = TRUE)
    df_selected <- df_merged[, c("V4", "MR1", "MR2","MR3", "MR61", "MR62", "MR63", "MR64", "MR65", "MR66", "MR67", "MR68", "MR69")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCounts.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Combine replicates (MR1, MR2, MR3) by averaging them
    df_merged$NT <- rowMeans(df_merged[, c("MR1", "MR2", "MR3")], na.rm = TRUE)
    df_merged$W8 <- rowMeans(df_merged[, c("MR61", "MR62", "MR63")], na.rm = TRUE)
    df_merged$W17 <- rowMeans(df_merged[, c("MR64", "MR65", "MR66")], na.rm = TRUE)
    df_merged$W28 <- rowMeans(df_merged[, c("MR67", "MR68", "MR69")], na.rm = TRUE)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # Log-transform the raw counts for each time point
    df_merged$NT <- log1p(df_merged$NT)
    df_merged$W8 <- log1p(df_merged$W8)
    df_merged$W17 <- log1p(df_merged$W17)
    df_merged$W28 <- log1p(df_merged$W28)
    df_selected <- df_merged[, c("V4", "NT", "W8","W17", "W28")]
    write.table(df_selected, paste0("./", folder, "/", folder, "_rawCountsAveraged_lgTransformed.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
}
