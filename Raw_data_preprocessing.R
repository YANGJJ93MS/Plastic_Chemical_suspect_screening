# Load necessary libraries
library(dplyr)
library(readr)


#####
  # Step 1: Load the data with all columns as characters to avoid parsing issues
  file_path <- "D:/UCSF_postdoc_topic/REVEAL_topics/MSDIAL_pipeline/Area_1_2024_11_03_20_04_54.txt" # Replace with your actual file path
  # data <- read_tsv(file_path, col_names = FALSE, skip_empty_rows = FALSE)
  data <- read_tsv(file_path, col_names = FALSE, col_types = cols(.default = "c"), skip_empty_rows = FALSE)
  
  # Step 2: Extract sample metadata
  sample_meta_raw <- data[1:4, ] # Selecting the first four rows
  # colnames(sample_meta) <- paste0("V", 1:ncol(sample_meta)) # Assign generic column names
  
  # Transpose and clean the metadata
  sample_meta <- sample_meta_raw %>%
    t() %>%                                      # Transpose the data to switch rows to columns
    as.data.frame() %>%                          # Convert to a data frame
    setNames(c("Class", "File type", "Injection order", "Batch ID")) %>% # Rename columns
    filter(!is.na(Class) & Class != "NA")                # Remove rows with all NA values
  sample_meta <- sample_meta[-1, ]  
  
  # Step 3: Extract peak area data
  peak_area_data <- data[-(1:4), ] # Removing the first four rows for peak area information
  colnames(peak_area_data) <- as.character(peak_area_data[1, ]) # Set column names from the first row
  peak_area_data <- peak_area_data[-1, ] # Remove the row used for column names



  
######
# import positive data
  file_path <- "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/Area_1_2024_11_04_23_46_26.txt"
  # data <- read_tsv(file_path, col_names = FALSE, skip_empty_rows = FALSE, show_col_types = FALSE)
  data <- read_tsv(file_path, col_names = FALSE, col_types = cols(.default = "c"), skip_empty_rows = FALSE)
  
  sample_meta_raw <- data[1:5, ] 
  # Transpose and clean the metadata
  sample_meta <- sample_meta_raw %>%
    t() %>%                                      # Transpose the data to switch rows to columns
    as.data.frame() %>%                          # Convert to a data frame
    setNames(c("Class", "File type", "Injection order", "Batch ID","sample_name")) %>% # Rename columns
    filter(!is.na(Class) & Class != "NA")                # Remove rows with all NA values
  sample_meta <- sample_meta[-1, ]  
  rownames(sample_meta) <- NULL
  
  # Step 3: Extract peak area data
  peak_area_data <- data[-(1:4), ] # Removing the first four rows for peak area information
  colnames(peak_area_data) <- as.character(peak_area_data[1, ]) # Set column names from the first row
  peak_area_data <- peak_area_data[-1, ] # Remove the row used for column names
  
  #drop the unwanted columns
  peak_area_data <- peak_area_data %>%
    select(`Average Rt(min)`, `Average Mz`, (which(names(peak_area_data) == "MS/MS spectrum")+1):(ncol(peak_area_data)-2))%>%
    mutate(across(everything(), as.numeric))

  #blank masking
  #remove features that are likely due to background ions and contaminants. intensity from samples 
  sample_id <- sample_meta %>%
    filter(`File type` == 'Sample')
  blank_id <- sample_meta %>%
    filter(`File type` == 'Blank') 
  
  #caculate the row means for smaples and blanks
  peak_area_with_means <- peak_area_data %>%
    mutate(
      sample_mean = rowMeans(select(., all_of(sample_id$sample_name)), na.rm=TRUE),
      blank_mean = rowMeans(select(., all_of(blank_id$sample_name)), na.rm=TRUE)
  )
  #filter rows where sample mean is at least 3 times the blank mean
  # filtered_peak_area <- peak_area_with_means %>%
  #   filter(sample_mean >= 3* blank_mean)
  filtered_peak_area <- peak_area_with_means %>%
    filter(sample_mean >= 0* blank_mean)
  metadata_cols <- filtered_peak_area %>%
    select(`Average Rt(min)`, `Average Mz`)
  
  #normalization based on internal standards maybe???
  #remove samples not needed for downstream analysis namely QC samples, blanks, etc.
  filtered_peak_sample <- filtered_peak_area %>%
    select(., all_of(sample_id$sample_name))
  
  #drop unfrequent columns, features; with high detection frequency
  filtered_peak_sample <- filtered_peak_sample %>%
    #detection frequency
    mutate(detection_frequency = rowSums(.!=0, na.rm=TRUE)/ncol(.)*100)%>%
    #fitlered features with >= 70%
    # filter(detection_frequency>=70)%>%
    filter(detection_frequency>=0)%>%
    select(-detection_frequency)
  
  result_df <- bind_cols(metadata_cols[rownames(filtered_peak_sample), ], filtered_peak_sample)
  sample_meta_filtered <- metadata_cols[rownames(filtered_peak_sample),]
  # write.csv(file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/filtered_peak_area_for_bc.csv", result_df)
  write.csv(file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/raw_peak_area_before_bc.csv", result_df)
  
  #impute column values, impute 0 value with row min
  
  # filtered_peak_sample <- filtered_peak_sample %>%
  #   rowwise()%>%
  #   mutate(across(everything(), ~ifelse(.==0, min(cur_data()[cur_data()!=0], na.rm=TRUE),.)))%>%
  #   ungroup()
  # 
  # filtered_results <- bind_cols(meta_data_cols, filtered_peak_sample)
  
  #use data.table
  library(data.table)
  setDT(filtered_peak_sample)
  #calculate row min of non-zeor vlaues
  row_mins <- apply(filtered_peak_sample, 1, function(row){
    min(row[row!=0], na.rm = TRUE)
  })
  
  #replace 0 with rowmin
  filtered_peak_sample <- as.data.table(
    t(apply(filtered_peak_sample,1,function(row,min_vals){
      row[row==0] <- min_vals
      return(row)
    }, min_vals = row_mins))
  )
  
  #Batch correction
  ##batch information
  ##peak area matrix
  batch <- sample_id$`Batch ID`
  count_matrix <- as.matrix(filtered_peak_sample)
  adjusted_counts <- ComBat_seq(count_matrix, batch = batch, group = NULL, full_mod = FALSE)
  save.image(file = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/pos_data_processing.RData' )
  
  ##plot PCA before batch correction and after batch correction
  library(ggplot2)
  pca_before <- prcomp(t(count_matrix), scale. = TRUE)
  
  pca_data_before <- data.frame(
    PC1=pca_before$x[,1],
    PC2 = pca_before$x[,2],
    Batch = factor(batch)
  )

  p_before <- ggplot(pca_data_before, aes(x =PC1, y=PC2, color = Batch))+
    geom_point(size=3)+
    labs(title='PCA Before Batch Correction', x ='PC1', y ='PC2', color = "Batch ID")+
    theme_minimal()
  p_before
  ggsave(p_before, file = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/PCA_before_batchcorr.tiff',
         dpi=300, height = 8, width=10)
  
  ##plot after batch correction
  pca_after <- prcomp(t(adjusted_counts), scale.=TRUE)
  pca_data_after <- data.frame(
    PC1 = pca_after$x[,1],
    PC2 = pca_after$x[,2],
    Batch = factor(batch)
  )
  p_after <- ggplot(pca_data_after, aes(x =PC1, y=PC2, color = Batch))+
    geom_point(size=3)+
    labs(title='PCA after Batch Correction', x ='PC1', y ='PC2',color = "Batch ID")+
    theme_minimal()
  ggsave(p_after, file = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/PCA_after_batchcorr.tiff',
         dpi=300, height = 8, width=10)
  
  #make ms1 search, and return the suspect list
  #combined meta column with adjusted counts
  result_df2 <- bind_cols(sample_meta_filtered, adjusted_counts)
  write.csv(file = "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/filtered_peak_area_after_bc.csv", result_df2)
  
  
  #Log transformation, correlation heatmap among features, or molecular network
  #
  
  
######
##prepare sample matrix for combat correction, demo code for combat_seq correction
##colname corresponds to each sample in each batch
##rows corresponds to the features 
##note, transpose to fit the gene expression format
##do not include condition for the correction this time
count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
batch <- c(rep(1, 4), rep(2, 4))
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL, full_mod=FALSE)


