##detection frequency for suspect features. 
suspect_identification <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/Final_Identification_pos_neg_v1.csv")
neg_peak_area <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Neg_AlignmentResults/filtered_peak_area_after_bc_10.csv")
pos_peak_area <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/filtered_peak_area_after_bc_10.csv")
demogra_data <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/demographic_data/Demoographics_fmt_copy.csv')
sample_name <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/demographic_data/sample_id_serum_clean.csv')

#combined the retention time and mz as the peak id
suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")

#add the education to the pptid
demogra_data <- demogra_data %>% 
  mutate(ga_weeks_mr = as.numeric(ga_weeks_mr),
         Terms = case_when(
           is.na(ga_weeks_mr) ~ NA_character_,
           ga_weeks_mr < 37 ~'preterm',
           ga_weeks_mr >= 37 ~'full_term'),
         GDM_diag = case_when(
           gstdiab_mr == "1, Yes" ~ "GDM",
           gstdiab_mr == "0, No"  ~ "Non-GDM",
           gstdiab_mr == "."      ~ NA_character_),
         Edu = case_when(
           edu_m_pre_cv == "5, Graduate Degree" ~ "Post_grad",
           edu_m_pre_cv %in% c("4, Bachelor's Degree", "3, Some College or AA Degree") ~ "College",
           edu_m_pre_cv %in% c("1, <High School", "2, HS Degree or GED") ~ "Highschool/less",
           edu_m_pre_cv %in% c(".", "") | is.na(edu_m_pre_cv) ~ NA_character_
         )) %>%
  select(ppt_id,Terms,GDM_diag, Edu, everything()) %>%
  mutate(ppt_id = as.character(ppt_id))


##plot the heat map of the each suspect identfication
##make log transformation
neg_peak_area_new <- neg_peak_area[,4:ncol(neg_peak_area)]
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
neg_peak_area_averaged$peak_id <- paste(neg_peak_area$`Average.Rt.min`, neg_peak_area$`Average.Mz`, sep = "_")

#log 2 transform for the columns value except 'peak_id'
#impute the na value with the minimum row value
cols_to_transform <- setdiff(names(neg_peak_area_averaged), "peak_id")

# 1. Log2 transform the participant columns
neg_peak_area_averaged[cols_to_transform] <- lapply(neg_peak_area_averaged[cols_to_transform], 
                                                    function(x) as.numeric(as.character(x)))
neg_peak_area_averaged[cols_to_transform] <- log2(neg_peak_area_averaged[cols_to_transform])

# 2. Impute NA values with the minimum row value for each row (excluding peak_id)
neg_peak_area_averaged[cols_to_transform] <- t(
  apply(neg_peak_area_averaged[cols_to_transform], 1, function(x) {
    if (all(is.na(x))) {
      # If an entire row is NA, just return it as is or decide on another strategy
      return(x)
    }
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
  })
)

#combined the retention time and mz as the peak id
suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")
neg_suspect <- suspect_identification%>%
  filter(polarity=='neg')

#subset the row according to the peak_id from suspect_identification data
neg_peak_area_averaged <- neg_peak_area_averaged%>%
  filter(peak_id %in% neg_suspect$peak_id)

##perform clustering plot








