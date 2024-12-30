#get detection frequency for 40 featuer in negative and 25 features in positive mode data
# import positive data
file_path <- "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/Area_1_2024_11_04_23_46_26.txt"
file_path2 <- 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Neg_AlignmentResults/Area_1_2024_11_04_16_24_07.txt'
suspect_identification <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/Combined_pos_neg_features_after_annotation_for_demographic.csv")
sample_name <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/demographic_data/sample_id_serum_clean.csv')

#combined the retention time and mz as the peak id
suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")

data <- read_tsv(file_path, col_names = FALSE, col_types = cols(.default="c"), skip_empty_rows = FALSE)
data2 <- read_tsv(file_path2, col_names = FALSE, col_types = cols(.default="c"), skip_empty_rows = FALSE)

sample_meta_raw <- data2[1:5, ] 
# Transpose and clean the metadata
sample_meta <- sample_meta_raw %>%
  t() %>%                                      # Transpose the data to switch rows to columns
  as.data.frame() %>%                          # Convert to a data frame
  setNames(c("Class", "File type", "Injection order", "Batch ID","sample_name")) %>% # Rename columns
  filter(!is.na(Class) & Class != "NA")                # Remove rows with all NA values
sample_meta <- sample_meta[-1, ]  
rownames(sample_meta) <- NULL

# Step 3: Extract peak area data
peak_area_data <- data2[-(1:4), ] # Removing the first four rows for peak area information
colnames(peak_area_data) <- as.character(peak_area_data[1, ]) # Set column names from the first row
peak_area_data <- peak_area_data[-1, ] # Remove the row used for column names
Adduct_list <- as.character(peak_area_data$`Adduct type`)
dim(peak_area_data)

#drop the unwanted columns
peak_area_data <- peak_area_data %>%
  select(`Average Rt(min)`, `Average Mz`,(which(names(peak_area_data) == "MS/MS spectrum")+1):(ncol(peak_area_data)-2))%>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(Adduct = Adduct_list) %>%
  # filter(Adduct%in%'[M+H]+') %>%
  filter(Adduct%in%'[M-H]-') %>%
  select(-c("Adduct"))
dim(peak_area_data)

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
filtered_peak_area <- peak_area_with_means %>%
  filter(sample_mean >= 3* blank_mean)

metadata_cols <- filtered_peak_area %>%
  select(`Average Rt(min)`, `Average Mz`)


#normalization based on internal standards maybe???
#remove samples not needed for downstream analysis namely QC samples, blanks, etc.
filtered_peak_sample <- filtered_peak_area %>%
  select(., all_of(sample_id$sample_name))

#drop unfrequent columns, features; with high detection frequency
filtered_peak_sample_0 <- filtered_peak_sample %>%
  #detection frequency
  mutate(detection_frequency = rowSums(.>5000, na.rm=TRUE)/ncol(.)*100)%>%
  filter(detection_frequency>=70)%>%
  select(-detection_frequency)


result_df_0 <- bind_cols(metadata_cols[rownames(filtered_peak_sample_0), ], filtered_peak_sample_0)

dim(result_df_0)

neg_peak_area_new <- result_df_0[,3:ncol(result_df_0)]
#average the peak area by duplicated column name, and assign the unique column name and the average peak area as the column value
#return the column
clean_names <- sub("BH([0-9]{5}).*", "\\1", colnames(neg_peak_area_new))

# Assign the cleaned names to the replicate dataframe
colnames(neg_peak_area_new) <- clean_names
neg_peak_area_new[] <- lapply(neg_peak_area_new, function(x) as.numeric(as.character(x)))

# Step 3: For each unique sample ID, average across all columns with that name
# Use sapply to loop over unique sample names and compute rowMeans
neg_peak_area_averaged <- as.data.frame(sapply(unique(clean_names), function(sample_id) {
  rowMeans(neg_peak_area_new[, colnames(neg_peak_area_new) == sample_id, drop = FALSE], na.rm = TRUE)
}))

colnames(neg_peak_area_averaged) <- unique(sample_name$RO1.REVEAL.Sample.number[!is.na(sample_name$RO1.REVEAL.Sample.number)])

neg_peak_area_averaged$peak_id <- paste(result_df_0$`Average Rt(min)`, result_df_0$`Average Mz`, sep = "_")


#combined the retention time and mz as the peak id
suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")

# neg_suspect <- suspect_identification%>%
#   filter(polarity=='pos')

neg_suspect <- suspect_identification%>%
  filter(polarity=='neg')
#subset the row according to the peak_id from suspect_identification data
neg_peak_area_df <- neg_peak_area_averaged%>%
  filter(peak_id %in% neg_suspect$peak_id)

#get the detected number of participant for each peak_id
#threshold 5000
count_over_5000 <- rowSums(neg_peak_area_df[,!(names(neg_peak_area_df)%in%"peak_id")] >0)
count_over_5000

neg_detectf <- data.frame(
  peak_id = neg_peak_area_df$peak_id,
  det_num_sample = count_over_5000,
  det_frequency = count_over_5000/200*100
)


pos_detectf <- data.frame(
  peak_id = neg_peak_area_df$peak_id,
  det_num_sample = count_over_5000,
  det_frequency = count_over_5000/200*100
)


new_combined = rbind(neg_detectf,pos_detectf)
dim(new_combined)
write.csv(new_combined, file = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/combined_detection_frequency.csv')
