library(dplyr)
library(tidyr)
library(ggplot2)
# install.packages('ggsignif')
library(ggsignif)

neg_peak_area <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Neg_AlignmentResults/filtered_peak_area_after_bc_70.csv")
pos_peak_area <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Pos_AlignmentResults/filtered_peak_area_after_bc_70.csv')
sample_name <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/demographic_data/sample_id_serum_clean.csv')
demogra_data <- read.csv('D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/demographic_data/Demoographics_fmt_copy.csv')
# suspect_identification <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/Final_Identification_pos_neg_v1.csv")
suspect_identification <- read.csv("D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/Combined_pos_neg_features_after_annotation_for_demographic.csv")


  #combined the retention time and mz as the peak id
  suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")

  #add preterm and term lable to the pptid
  #defined the preterm group, term group, GDM group, DM group, and the education group
  #add the GDM label to the pptid
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
             edu_m_pre_cv %in% c(".", "") | is.na(edu_m_pre_cv) ~ NA_character_),
           Race = case_when(
             race_m_cv == '1, White'~'NH-White',
             race_m_cv == '2, Black or African American' ~ 'NH-Black',
             race_m_cv == '3, Asian' ~ 'NH-Asian',
             race_m_cv %in% c('4, Native Hawaiian or Other Pacific Islander',
                           '5, American Indian or Alaska Native','6, More than 1 race',
                           '-8, Unknown, Not Reported or Other') ~ 'other',
             race_m_cv == 'hispanic' ~ 'Latina'
           ),
           Nativity = case_when(
             usborn_m_cv == '1, US Born' ~ 'US',
             usborn_m_cv %in% c('.','0, Foreign Born') ~'other'
           ),
           Smoke = case_when(
             smoke_cv =='2, Former smoker' ~'smoke',
             smoke_cv == '1, Never smoked' ~'non-smoke',
             smoke_cv == '.' ~ NA_character_
           )
           ) %>%
    select(ppt_id,Terms,GDM_diag, Edu, everything())
  
  demogra_data <- demogra_data %>%
    mutate(ppt_id = as.character(ppt_id))
  
  print(dim(demogra_data[demogra_data$Terms=='preterm',]))
  print(dim(demogra_data[demogra_data$Terms=='full_term',]))
  dim(demogra_data)
  
  
  #######################
  ##handling negative peak area
  #reformat the peak area data
  #match the participant id with the sample name in peak area
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
  # Make sure these columns are numeric
  neg_peak_area_averaged[cols_to_transform] <- lapply(neg_peak_area_averaged[cols_to_transform], 
                                                      function(x) as.numeric(as.character(x)))
  
  neg_peak_area_averaged[cols_to_transform] <- log2(neg_peak_area_averaged[cols_to_transform])
  
  # 2. Impute NA values with the minimum row value for each row (excluding peak_id)
  neg_peak_area_averaged[cols_to_transform] <- t(apply(neg_peak_area_averaged[cols_to_transform], 1, function(x){
      if (all(is.na(x))) {
        # If an entire row is NA
        return(x)
      }
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    }))
  

  #combined the retention time and mz as the peak id
  suspect_identification$peak_id <- paste(suspect_identification$`Average.Rt.min`, suspect_identification$`Average.Mz`, sep = "_")
  neg_suspect <- suspect_identification%>%
    filter(polarity=='neg')
  
  #subset the row according to the peak_id from suspect_identification data
  neg_peak_area_averaged <- neg_peak_area_averaged%>%
    filter(peak_id %in% neg_suspect$peak_id)
  
 
  ##################################################
  ##select column without NA value accroding to the demogra_data when perform test for each variable, ppt_id~terms, GDM_diag, Edu
  ##perform Krushkal-wallis rank sum test to test if the import feature that are highly impact by the label.
  ## set the p-value threshold to 0.05.
  #make the box plot of the relative concentration for the peak_id with signification p-value < 0.05
  # Convert peak_area_data to a long format and merge with demographic data
  # Ensure ppt_id is character type
  
  # Reshape peak area data
  peak_area_long <- neg_peak_area_averaged %>%
    pivot_longer(-peak_id, names_to = "ppt_id", values_to = "peak_area") 
  
  peak_area_long <- peak_area_long%>%
    left_join(demogra_data, by = "ppt_id")
  
  
  # Step 2: Perform Kruskal-Wallis tests with multiple comparison adjustment
  results <- list()  # To store results
  # demographic_vars <- c("Terms", "GDM_diag", "Edu")  # Variables to test
  demographic_vars <- c("Terms", "GDM_diag", "Edu",'Race','Nativity')
  # demographic_vars <- c("Terms",'GDM_diag')
  
  # Perform Kruskal-Wallis test and adjust p-values
  for (var in demographic_vars) {
    # Remove rows with NA in the demographic variable
    filtered_data <- peak_area_long %>%
      filter(!is.na(.data[[var]]))
    
    # Perform Kruskal-Wallis test for each peak_id
    test_results <- filtered_data %>%
      group_by(peak_id) %>%
      summarise(
        p_value = kruskal.test(peak_area ~ .data[[var]])$p.value
      ) %>%
      mutate(
        variable = var,
        adj_p_value = p.adjust(p_value, method = "BH") # Adjust p-values
      )
    
    results[[var]] <- test_results
  }
  
  # Combine results for all variables
  all_results <- bind_rows(results)
  
  # Step 3: Filter for significant adjusted p-values
  significant_results <- all_results %>%
    filter(p_value < 0.05)
  message("Significant results with p-value < 0.05: ", nrow(significant_results))
  
  write.csv(significant_results, file='D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_negative/significant_results_p0.05_20241229.csv')
  
  # Step 4: Generate Box Plots with P-values
  filepath <- 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_negative/'
  # dir.create(paste0(filepath,'negative_112chemicals_20241229/'))
  filepath <- paste0(filepath,'negative_112chemicals_20241229/')
  significant_peak_ids <- significant_results$peak_id
  
  for (var in demographic_vars) {
    for (peak_id in significant_peak_ids) {
      plot_data <- peak_area_long %>%
        filter(peak_id == !!peak_id, !is.na(.data[[var]]))
      
      # Extract p-values for the current peak_id and variable
      p_info <- significant_results %>%
        filter(peak_id == !!peak_id, variable == !!var) %>%
        select(p_value, adj_p_value)
      
      if (nrow(p_info) == 0) next
      
      raw_p <- formatC(p_info$p_value, format = "e", digits = 3)
      adj_p <- formatC(p_info$adj_p_value, format = "e", digits = 3)
      
      if (var =='Terms')
        comparisons = list(c("full_term", "preterm"))
      else if (var == 'GDM_diag')
        comparisons = list(c('Non-GDM','GDM'))
      else if (var == 'Edu')
        comparisons =  list(
          c("College", "Highschool/less"),
          c("College", "Post_grad"),
          c("Highschool/less", "Post_grad")
        )
      else if (var =='Nativity')
        comparisons = list(c('US','other'))
      
      else if (var=='Race')
        comparisons = list(
          c('other', 'Latina'),
          c('other', 'Latina'),
          c('other', 'NH-Black'),
          c('other','NH-White'),
          c('other','NH-Asian'),
          c('Latina','NH-Black'),
          c('Latina','NH-White'),
          c('Latina','NH-Asian'),
          c('NH-White','NH-Black'),
          c('NH-White','NH-Asian'),
          c('NH-Black','NH-Asian')
        )
        
      
      # Generate the box plot
      p <- ggplot(plot_data, aes(x = .data[[var]], y = peak_area)) +
        
        # Add boxplot with detailed explanation
        geom_boxplot(outlier.shape = 21, outlier.color = "red", outlier.fill = "white", 
                     outlier.size = 4, size = 1.5,
                     coef = 1.5) +
        
        # Add jitter points to show individual data points
        # geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
        geom_signif(
          comparisons = comparisons,  # Specify the conditions to compare
          map_signif_level = TRUE,  # Automatically map p-values to asterisks
          step_increase = 0.05,  # Adjust vertical position of the significance bar
          y_position =max(plot_data$peak_area, na.rm = TRUE) + 0.2* max(plot_data$peak_area, na.rm = TRUE),
          textsize = 4,  # Adjust size of asterisk
          size = 1 # adjust the sig lines thickness
        )+
        
        # Customize plot labels
        labs(
          title = paste("Peak:", peak_id, "| p:", raw_p),  # Display peak ID and p-values in the title
          x = var,
          y = "Relative Concentration (Peak Area)"
        ) +
        
        # Customize the theme for bold fonts and larger numeric labels
        theme_minimal(base_size = 18) +  # Set base font size
        theme(
          # Title
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),  # Large bold title
          
          # Axis Titles
          axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # Large bold x-axis label
          axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Large bold y-axis label
          
          # Axis Text (Numeric Labels)
          axis.text.x = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold x-axis text
          axis.text.y = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold y-axis text
          
          # Legend
          legend.title = element_text(size = 18, face = "bold", color = "black"),  # Bold legend title
          legend.text = element_text(size = 16, color = "black"),  # Legend text
          
          # Border and Grid
          panel.border = element_rect(color = "black", fill = NA, size = 3),  # Add black border
          panel.grid.major = element_line(color = "gray90", size = 1),  # Subtle grid lines
          panel.grid.minor = element_blank()  # Remove minor grid lines
        )
      
      
      # Save the plot in .tif format
      ggsave(
        filename = paste0(filepath, "boxplot_sig", peak_id, "_", var, ".tif"),
        plot = p, width = 8, height = 10, dpi = 300, device = "tiff"
      )
      
      message("Saved plot for Peak ID:", peak_id, " | Variable:", var)
    }
  }
  
  
  ########################################
  library(data.table)
  neg_peak_area_pca <- neg_peak_area_averaged%>%
    filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='Terms']))
  # neg_peak_area_pca <- neg_peak_area_averaged%>%
  #   filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='GDM_diag']))
  dt <- as.data.table(neg_peak_area_pca)
  pptid <- colnames(dt)[1:200]
  peakid <- dt$peak_id
  tran_dat <- setDT(data.table::transpose(dt[,-which(names(dt)=='peak_id'), with=FALSE]))
  tran_dat$ppt_id <- pptid
  tran_dat <- merge(x=tran_dat, y = demogra_data[, c('ppt_id','Terms')], by = 'ppt_id', all.x = TRUE)
  # tran_dat <- merge(x=tran_dat, y = demogra_data[, c('ppt_id','GDM_diag')], by = 'ppt_id', all.x = TRUE)
  # tran_dat <- tran_dat[!is.na(GDM_diag)]
  tran_dat <- tran_dat[!is.na(Terms)]
  pca_dat <- as.data.frame(scale(tran_dat[,-c('ppt_id','Terms')], center=TRUE, scale=TRUE))
  # pca_dat <- as.data.frame(scale(tran_dat[,-c('ppt_id','Edu')], center=TRUE, scale=TRUE))
  # pca_dat <- as.data.frame(scale(tran_dat[,-c('ppt_id','GDM_diag')], center=TRUE, scale=TRUE))
  rownames(pca_dat) <- tran_dat$ppt_id
  pca_res <- prcomp(pca_dat, scale. = TRUE)
  pc1_scores <-as.data.frame(pca_res$x)
  # pc1_scores$GDM_diag <- tran_dat$GDM_diag
  pc1_scores$Terms <- tran_dat$Terms
  # ttest_res <- t.test(PC1 ~ GDM_diag, data = pc1_scores)
  ttest_res <- t.test(PC1 ~ Terms, data = pc1_scores)
  p_value <- ttest_res$p.value
  p_value
  # Define significance level based on p-value
  significance <- ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01, "**",
                                ifelse(p_value < 0.05, "*", "ns")))
  
  # Add box plot with p-value and significance level
  p <- ggplot(pc1_scores, aes(x = Terms, y = PC1, fill = Terms)) +
    geom_boxplot(outlier.shape = 21, outlier.color = "red", outlier.fill = "white", 
                 outlier.size = 4, size = 1.5,
                 coef = 1.5) +
    labs(
      title = paste0("PC1 Scores by","Terms"),
      x = "Terms",
      y = "PC1 Scores"
    ) +
    # Add significance line using ggsignif
    geom_signif(
      # comparisons = list(c("GDM", "Non-GDM")),  # Specify the groups to compare
      comparisons = list(c("preterm", "full_term")),  # Specify the groups to compare
      annotations = significance,                    # Add the significance level
      map_signif_level = TRUE,                       # Automatically maps p-value
      y_position = max(pc1_scores$PC1, na.rm = TRUE) * 1.05,  # Position above the boxes
      textsize = 6                                   # Size of significance annotation
    ) +
    # theme_minimal(base_size = 14) +
    # theme(
    #   plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    #   axis.title.x = element_text(size = 18, face = "bold"),
    #   axis.title.y = element_text(size = 18, face = "bold"),
    #   axis.text.x = element_text(size = 16, face = "bold"),
    #   axis.text.y = element_text(size = 16, face = "bold"),
    #   legend.position = "none"
    # )
    theme_minimal(base_size = 14) +  # Set base font size
    theme(
      # Title
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),  # Large bold title
      
      # Axis Titles
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # Large bold x-axis label
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Large bold y-axis label
      
      # Axis Text (Numeric Labels)
      axis.text.x = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold x-axis text
      axis.text.y = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold y-axis text
      
      # Legend
      legend.title = element_text(size = 18, face = "bold", color = "black"),  # Bold legend title
      legend.text = element_text(size = 16, color = "black"),  # Legend text
      
      # Border and Grid
      panel.border = element_rect(color = "black", fill = NA, size = 3),  # Add black border
      panel.grid.major = element_line(color = "gray90", size = 1),  # Subtle grid lines
      panel.grid.minor = element_blank()  # Remove minor grid lines
    )
  p
  # ggsave(filename = paste0(filepath, "PCA_loading_Plot_Term_112chem", ".tif"), plot = p, width = 8, height = 10, dpi = 300)
  ggsave(filename = paste0(filepath, "PCA_loading_Plot_Term_112chem", ".tif"), plot = p, width = 8, height = 10, dpi = 300)
  
  
  ###############################################
  #perform clustering with the deographic data and the peak area for each participant for each peak-id
  #make the clustering plot
  # Clustering Analysis Function
  library(cluster)
  library(pheatmap)
  library(tidyr)
  
  ##########
  # Step 1: Filter out participants with missing demographic data
  # Step 1: Prepare filtered demographic data
  # filtered_demogra_data <- demogra_data %>%
  #   mutate(ppt_id = as.character(ppt_id)) %>%
  #   filter(!is.na(Terms) & !is.na(GDM_diag) & !is.na(Edu))
  # filtered_demogra_data <- demogra_data %>%
  #   mutate(ppt_id = as.character(ppt_id))
  # 
  # # Filter peak area data to include only participants with complete demographic data
  # filtered_peak_area_long <- neg_peak_area_averaged %>%
  #   pivot_longer(-peak_id, names_to = "ppt_id", values_to = "peak_area") %>%
  #   left_join(filtered_demogra_data, by = "ppt_id")

  
  # Filepath for saving plots
  output_dir <- "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_negative/negative_112chemicals_20241229/"
  
  # Step 2: Define demographic variables for analysis
  # demographic_vars <- c("Terms")
  # demographic_vars <- c("GDM_diag") 
  demographic_vars <- c("Terms")
  
  # filtered_demogra_data <- demogra_data %>%
  #   mutate(ppt_id = as.character(ppt_id))
  # 
  # # Filter peak area data to include only participants with complete demographic data
  # filtered_peak_area_long <- pos_peak_area_averaged %>%
  #   pivot_longer(-peak_id, names_to = "ppt_id", values_to = "peak_area") %>%
  #   left_join(filtered_demogra_data, by = "ppt_id")
  
  #select features based on the sig results
  peak_area_long2 <- peak_area_long%>%
    filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='Terms']))
  
  # peak_area_long3 <- peak_area_long%>%
  #   filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='GDM_diag']))
  # 
  # peak_area_long4 <- peak_area_long%>%
  #   filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='Edu']))
  
  peak_area_long2 <- peak_area_long2 %>%
    filter(!is.na(.data[[demographic_vars]])) %>%
    select(peak_id,ppt_id,peak_area,Terms)
  # 
  # data_wide <- peak_area_long2 %>%
  #   pivot_wider(names_from = ppt_id, values_from = peak_area) %>%
  #   column_to_rownames("peak_id")
  # 
  # scaled_data <- scale(data_wide, center = TRUE, scale = TRUE)
  # 
  # scaled_data2 <- data_wide %>%
  #   rowwise() %>%
  #   mutate(across(everything(), ~ (. - mean(c_across(everything()))) / sd(c_across(everything())))) %>%
  #   ungroup()
    
  
  # Step 3: Hierarchical Clustering Function 
  perform_hierarchical_clustering <- function(data, demo_var) {
    
    # Filter data for the specific demographic variable
    data_filtered <- data %>%
      filter(!is.na(.data[[demo_var]]))
    
    # Transform data into wide format for clustering
    data_wide <- data_filtered %>%
      select(peak_id, ppt_id, peak_area)%>%
      pivot_wider(names_from = ppt_id, values_from = peak_area) %>%
      column_to_rownames("peak_id")
    
    # Scale the data, per aprticipant
    # scaled_data <- scale(data_wide, center = TRUE, scale = TRUE)
    
    scaled_data <-data_wide %>%
      as.data.frame()%>%
      rowwise() %>%
      mutate(across(everything(), ~ (. - mean(c_across(everything()))) / sd(c_across(everything())))) %>%
      ungroup() %>%
      as.data.frame()
    
    rownames(scaled_data) <-rownames(data_wide)
    
    # Perform Hierarchical Clustering
    dist_matrix <- dist(scaled_data, method = "euclidean")  # Compute distance matrix
    hc <- hclust(dist_matrix, method = "complete")          # Perform hierarchical clustering
    
    # Cut the den drogram into clusters
    clusters <- cutree(hc, k = 2)  # Cut into 4 clusters
    
    # Add cluster assignments to the data
    cluster_assignments <- data.frame(
      peak_id = rownames(scaled_data),
      cluster = as.factor(clusters)
    )
    
    # Prepare annotations using ppt_id
    demo_annotations <- data_filtered %>%
      distinct(ppt_id, Terms) %>%
      column_to_rownames("ppt_id")
    # demo_var <- demo_var
    # demo_annotations <- data_filtered %>%
    #   distinct(ppt_id, GDM_diag) %>%
    #   column_to_rownames("ppt_id")
    # demo_annotations <- data_filtered %>%
    #   distinct(ppt_id, Edu) %>%
    #   column_to_rownames("ppt_id")
    
    # Reorder columns based on group order
    # column_order <- demo_annotations %>%
    #   arrange(Edu) %>%
    #   rownames()
    # scaled_data <- scaled_data[, column_order]
    # Heatmap Visualization
    heatmap_file <- paste0(output_dir, "Heatmap_", demo_var, ".tif")
    png(heatmap_file, width = 1200, height = 1000, res = 200)
    pheatmap(
      scaled_data,
      cluster_rows = TRUE,            # Use hierarchical clustering for rows
      cluster_cols = TRUE,  # Use hierarchical clustering for columns
      # cluster_cols = FALSE,
      annotation_col = demo_annotations, # Add demographic annotations
      main = paste("Clustering Heatmap -", demo_var),
      color = colorRampPalette(c("blue", "white", "red"))(50),
      fontsize_row = 12,
      show_colnames = FALSE  # Hide participant IDs at the bottom
    )
    dev.off()
    
    
    ## Perform PCA 
    pca_result <- prcomp(scaled_data, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x) %>%
      rownames_to_column("peak_id") %>%
      left_join(cluster_assignments, by = "peak_id")

    ##perform box plot using PC 1 loadings
    # PCA Plot (without participant ID labels)
    pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
      geom_point(size = 4, alpha = 0.7) +
      labs(
        title = paste("PCA Clustering by", demo_var),
        x = "PC1 (Principal Component 1)",
        y = "PC2 (Principal Component 2)",
        color = "Cluster"
      ) +
      scale_color_manual(values = c("red", "blue", "green", "purple")) +
      theme_minimal()
    #
    # # Save the PCA plot
    ggsave(
      filename = paste0(output_dir, "PCA_Plot_", demo_var, ".tif"),
      plot = pca_plot,
      width = 8,
      height = 6,
      dpi = 300
    )

    # print(pca_plot)

    
    return(list(
      pca_plot = pca_plot,
      hc_result = hc,
      cluster_assignments = cluster_assignments
    ))
  }
  
  
  # step: Hierarchical clustering by compounds
  
  
  # Step 4: Run hierarchical clustering for each demographic variable
  results <- list()
  for (var in demographic_vars) {
    message("Running hierarchical clustering for:", var)
    results[[var]] <- perform_hierarchical_clustering(peak_area_long2, var)
  }
  
  
  
  #########################################################################################################
  ##handling positive data
  pos_peak_area_new <- pos_peak_area[,4:ncol(pos_peak_area)]
  #log 2 transform for the columns value except 'peak_id'
  #impute the na value with the minimum row value
  #average the peak area by duplicated column name, and assign the unique column name and the average peak area as the column value
  clean_names <- sub("BH([0-9]{5}).*", "\\1", colnames(pos_peak_area_new))
  
  # Assign the cleaned names to the replicate dataframe
  colnames(pos_peak_area_new) <- clean_names
  pos_peak_area_new[] <- lapply(pos_peak_area_new, function(x) as.numeric(as.character(x)))
  
  # Step 3: For each unique sample ID, average across all columns with that name
  # Use sapply to loop over unique sample names and compute rowMeans
  pos_peak_area_averaged <- as.data.frame(sapply(unique(clean_names), function(sample_id) {
    rowMeans(pos_peak_area_new[, colnames(pos_peak_area_new) == sample_id, drop = FALSE], na.rm = TRUE)
  }))
  
  colnames(pos_peak_area_averaged) <- unique(sample_name$RO1.REVEAL.Sample.number[!is.na(sample_name$RO1.REVEAL.Sample.number)])
  pos_peak_area_averaged$peak_id <- paste(pos_peak_area$`Average.Rt.min`, pos_peak_area$`Average.Mz`, sep = "_")
  
  #log 2 transform for the columns value except 'peak_id'
  #impute the na value with the minimum row value
  cols_to_transform <- setdiff(names(pos_peak_area_averaged), "peak_id")
  
  # 1. Log2 transform the participant columns
  # Make sure these columns are numeric
  pos_peak_area_averaged[cols_to_transform] <- lapply(pos_peak_area_averaged[cols_to_transform], 
                                                      function(x) as.numeric(as.character(x)))
  
  pos_peak_area_averaged[cols_to_transform] <- log2(pos_peak_area_averaged[cols_to_transform])
  
  # 2. Impute NA values with the minimum row value for each row (excluding peak_id)
  pos_peak_area_averaged[cols_to_transform] <- t(
    apply(pos_peak_area_averaged[cols_to_transform], 1, function(x) {
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
  pos_suspect <- suspect_identification%>%
    filter(polarity=='pos')
  
  #subset the row according to the peak_id from suspect_identification data
  pos_peak_area_averaged <- pos_peak_area_averaged%>%
    filter(peak_id %in% pos_suspect$peak_id)
  peak_data <- pos_peak_area_averaged %>%
    select(-peak_id)
  
  # Reshape peak area data
  peak_area_long <- pos_peak_area_averaged %>%
    pivot_longer(-peak_id, names_to = "ppt_id", values_to = "peak_area") 
  peak_area_long <- peak_area_long%>%
    left_join(demogra_data, by = "ppt_id")
  
  # Step 2: Perform Kruskal-Wallis tests with multiple comparison adjustment
  results <- list()  # To store results
  # demographic_vars <- c("Terms", "GDM_diag")  # Variables to test
  demographic_vars <- c("Terms", "GDM_diag", "Edu",'Race','Nativity')
  
  # Perform Kruskal-Wallis test and adjust p-values
  for (var in demographic_vars) {
    # Remove rows with NA in the demographic variable
    filtered_data <- peak_area_long %>%
      filter(!is.na(.data[[var]]))
    
    # Perform Kruskal-Wallis test for each peak_id
    test_results <- filtered_data %>%
      group_by(peak_id) %>%
      summarise(
        p_value = kruskal.test(peak_area ~ .data[[var]])$p.value
      ) %>%
      mutate(
        variable = var,
        adj_p_value = p.adjust(p_value, method = "BH") # Adjust p-values
      )
    
    results[[var]] <- test_results
  }
  
  # Combine results for all variables
  all_results <- bind_rows(results)
  
  # Step 3: Filter for significant adjusted p-values
  significant_results <- all_results %>%
    filter(p_value < 0.05)
  
  message("Significant results with p-value < 0.05: ", nrow(significant_results))
  write.csv(significant_results, file='D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_positive/112chemical_significant_results_p0.05_20291229.csv')
  
  
  # Step 4: Generate Box Plots with P-values
  filepath <- 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_positive/'
  dir.create(paste0(filepath,'positive_112chemical_20241229/'))
  filepath <- paste0(filepath,'positive_112chemicals_20241229/')
  significant_peak_ids <- significant_results$peak_id
  
  # for (var in demographic_vars) {
  #   for (peak_id in significant_peak_ids) {
  #     plot_data <- peak_area_long %>%
  #       filter(peak_id == !!peak_id, !is.na(.data[[var]]))
  #     
  #     # Extract p-values for the current peak_id and variable
  #     p_info <- significant_results %>%
  #       filter(peak_id == !!peak_id, variable == !!var) %>%
  #       select(p_value, adj_p_value)
  #     
  #     if (nrow(p_info) == 0) next
  #     
  #     raw_p <- formatC(p_info$p_value, format = "e", digits = 3)
  #     adj_p <- formatC(p_info$adj_p_value, format = "e", digits = 3)
  #     
  #     if (var =='Terms')
  #       comparisons = list(c("full_term", "preterm"))
  #     else if (var == 'GDM_diag')
  #       comparisons = list(c('Non-GDM','GDM'))
  #     
  #     # Generate the box plot
  #     p <- ggplot(plot_data, aes(x = .data[[var]], y = peak_area)) +
  #       
  #       # Add boxplot with detailed explanation
  #       geom_boxplot(outlier.shape = 21, outlier.color = "red", outlier.fill = "white", 
  #                    outlier.size = 4, size = 1.5,
  #                    coef = 1.5) +
  #       
  #       # Add jitter points to show individual data points
  #       # geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  #       geom_signif(
  #         comparisons = comparisons,  # Specify the conditions to compare
  #         map_signif_level = TRUE,  # Automatically map p-values to asterisks
  #         step_increase = 0.05,  # Adjust vertical position of the significance bar
  #         y_position =max(plot_data$peak_area, na.rm = TRUE) + 0.2* max(plot_data$peak_area, na.rm = TRUE),
  #         textsize = 12,  # Adjust size of asterisk
  #         size = 1 # adjust the sig lines thickness
  #       )+
  #       
  #       # Customize plot labels
  #       labs(
  #         title = paste("Peak:", peak_id, "| p:", raw_p),  # Display peak ID and p-values in the title
  #         x = var,
  #         y = "Relative Concentration (Peak Area)"
  #       ) +
  #       
  #       # Customize the theme for bold fonts and larger numeric labels
  #       theme_minimal(base_size = 18) +  # Set base font size
  #       theme(
  #         # Title
  #         plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),  # Large bold title
  #         
  #         # Axis Titles
  #         axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # Large bold x-axis label
  #         axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Large bold y-axis label
  #         
  #         # Axis Text (Numeric Labels)
  #         axis.text.x = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold x-axis text
  #         axis.text.y = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold y-axis text
  #         
  #         # Legend
  #         legend.title = element_text(size = 18, face = "bold", color = "black"),  # Bold legend title
  #         legend.text = element_text(size = 16, color = "black"),  # Legend text
  #         
  #         # Border and Grid
  #         panel.border = element_rect(color = "black", fill = NA, size = 3),  # Add black border
  #         panel.grid.major = element_line(color = "gray90", size = 1),  # Subtle grid lines
  #         panel.grid.minor = element_blank()  # Remove minor grid lines
  #       )
  #     
  #     
  #     # Save the plot in .tif format
  #     ggsave(
  #       filename = paste0(filepath, "boxplot_sig", peak_id, "_", var, ".tif"),
  #       plot = p, width = 8, height = 10, dpi = 300, device = "tiff"
  #     )
  #     
  #     message("Saved plot for Peak ID:", peak_id, " | Variable:", var)
  #   }
  # }
  
  for (var in demographic_vars) {
    for (peak_id in significant_peak_ids) {
      plot_data <- peak_area_long %>%
        filter(peak_id == !!peak_id, !is.na(.data[[var]]))
      
      # Extract p-values for the current peak_id and variable
      p_info <- significant_results %>%
        filter(peak_id == !!peak_id, variable == !!var) %>%
        select(p_value, adj_p_value)
      
      if (nrow(p_info) == 0) next
      
      raw_p <- formatC(p_info$p_value, format = "e", digits = 3)
      adj_p <- formatC(p_info$adj_p_value, format = "e", digits = 3)
      
      if (var =='Terms')
        comparisons = list(c("full_term", "preterm"))
      else if (var == 'GDM_diag')
        comparisons = list(c('Non-GDM','GDM'))
      else if (var == 'Edu')
        comparisons =  list(
          c("College", "Highschool/less"),
          c("College", "Post_grad"),
          c("Highschool/less", "Post_grad")
        )
      else if (var =='Nativity')
        comparisons = list(c('US','other'))
      
      else if (var=='Race')
        comparisons = list(
          c('other', 'Latina'),
          c('other', 'Latina'),
          c('other', 'NH-Black'),
          c('other','NH-White'),
          c('other','NH-Asian'),
          c('Latina','NH-Black'),
          c('Latina','NH-White'),
          c('Latina','NH-Asian'),
          c('NH-White','NH-Black'),
          c('NH-White','NH-Asian'),
          c('NH-Black','NH-Asian')
        )
      
      
      # Generate the box plot
      p <- ggplot(plot_data, aes(x = .data[[var]], y = peak_area)) +
        
        # Add boxplot with detailed explanation
        geom_boxplot(outlier.shape = 21, outlier.color = "red", outlier.fill = "white", 
                     outlier.size = 4, size = 1.5,
                     coef = 1.5) +
        
        # Add jitter points to show individual data points
        # geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
        geom_signif(
          comparisons = comparisons,  # Specify the conditions to compare
          map_signif_level = TRUE,  # Automatically map p-values to asterisks
          step_increase = 0.05,  # Adjust vertical position of the significance bar
          y_position =max(plot_data$peak_area, na.rm = TRUE) + 0.2* max(plot_data$peak_area, na.rm = TRUE),
          textsize = 4,  # Adjust size of asterisk
          size = 1 # adjust the sig lines thickness
        )+
        
        # Customize plot labels
        labs(
          title = paste("Peak:", peak_id, "| p:", raw_p),  # Display peak ID and p-values in the title
          x = var,
          y = "Relative Concentration (Peak Area)"
        ) +
        
        # Customize the theme for bold fonts and larger numeric labels
        theme_minimal(base_size = 18) +  # Set base font size
        theme(
          # Title
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),  # Large bold title
          
          # Axis Titles
          axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # Large bold x-axis label
          axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Large bold y-axis label
          
          # Axis Text (Numeric Labels)
          axis.text.x = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold x-axis text
          axis.text.y = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold y-axis text
          
          # Legend
          legend.title = element_text(size = 18, face = "bold", color = "black"),  # Bold legend title
          legend.text = element_text(size = 16, color = "black"),  # Legend text
          
          # Border and Grid
          panel.border = element_rect(color = "black", fill = NA, size = 3),  # Add black border
          panel.grid.major = element_line(color = "gray90", size = 1),  # Subtle grid lines
          panel.grid.minor = element_blank()  # Remove minor grid lines
        )
      
      
      # Save the plot in .tif format
      ggsave(
        filename = paste0(filepath, "boxplot_sig", peak_id, "_", var, ".tif"),
        plot = p, width = 8, height = 10, dpi = 300, device = "tiff"
      )
      
      message("Saved plot for Peak ID:", peak_id, " | Variable:", var)
    }
  }
  
  #########################
  library(data.table)
  pos_peak_area_averaged <- pos_peak_area_averaged%>%
    filter(peak_id%in%(significant_results$peak_id[significant_results$variable=='Terms']))
  dt <- as.data.table(pos_peak_area_averaged)
  pptid <- colnames(dt)[1:200]
  peakid <- dt$peak_id
  tran_dat <- setDT(data.table::transpose(dt[,-which(names(dt)=='peak_id'), with=FALSE]))
  tran_dat$ppt_id <- pptid
  tran_dat <- merge(x=tran_dat, y = demogra_data[, c('ppt_id','Terms')], by = 'ppt_id', all.x = TRUE)
  tran_dat <- tran_dat[!is.na(Terms)]
  pca_dat <- as.data.frame(scale(tran_dat[,-c('ppt_id','Terms')], center=TRUE, scale=TRUE))
  rownames(pca_dat) <- tran_dat$ppt_id
  pca_res <- prcomp(pca_dat, scale. = TRUE)
  pc1_scores <-as.data.frame(pca_res$x)
  pc1_scores$Terms <- tran_dat$Terms
  ttest_res <- t.test(PC1 ~ Terms, data = pc1_scores)
  p_value <- ttest_res$p.value
  p_value
  # Define significance level based on p-value
  significance <- ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01, "**",
                                ifelse(p_value < 0.05, "*", "ns")))
  
  # Add box plot with p-value and significance level
  p <- ggplot(pc1_scores, aes(x = Terms, y = PC1, fill = Terms)) +
    geom_boxplot(outlier.shape = 21, outlier.color = "red", outlier.fill = "white", 
                 outlier.size = 4, size = 1.5,
                 coef = 1.5) +
    labs(
      title = "PC1 Scores by Term",
      x = "Term",
      y = "PC1 Scores"
    ) +
    # Add significance line using ggsignif
    geom_signif(
      comparisons = list(c("full_term", "preterm")),  # Specify the groups to compare
      annotations = significance,                    # Add the significance level
      map_signif_level = TRUE,                       # Automatically maps p-value
      y_position = max(pc1_scores$PC1, na.rm = TRUE) * 1.05,  # Position above the boxes
      textsize = 6                                   # Size of significance annotation
    ) +
    # theme_minimal(base_size = 14) +
    # theme(
    #   plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    #   axis.title.x = element_text(size = 18, face = "bold"),
    #   axis.title.y = element_text(size = 18, face = "bold"),
    #   axis.text.x = element_text(size = 16, face = "bold"),
    #   axis.text.y = element_text(size = 16, face = "bold"),
    #   legend.position = "none"
    # )
    theme_minimal(base_size = 14) +  # Set base font size
    theme(
      # Title
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),  # Large bold title
      
      # Axis Titles
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # Large bold x-axis label
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Large bold y-axis label
      
      # Axis Text (Numeric Labels)
      axis.text.x = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold x-axis text
      axis.text.y = element_text(size = 18, face = "bold", color = "black"),  # Enlarged bold y-axis text
      
      # Legend
      legend.title = element_text(size = 18, face = "bold", color = "black"),  # Bold legend title
      legend.text = element_text(size = 16, color = "black"),  # Legend text
      
      # Border and Grid
      panel.border = element_rect(color = "black", fill = NA, size = 3),  # Add black border
      panel.grid.major = element_line(color = "gray90", size = 1),  # Subtle grid lines
      panel.grid.minor = element_blank()  # Remove minor grid lines
    )
  p
  filepath = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_positive/positive_112chemicals_20241229/'
  ggsave(filename = paste0(filepath, "PCA_loading_Plot_Terms_112chem_20241229", ".tif"), plot = p, width = 8, height = 10, dpi = 300)
  #########################
  
  
  library(tibble)
  library(pheatmap)
  
  ##other plot:
  # filtered_demogra_data <- demogra_data %>%
  #   mutate(ppt_id = as.character(ppt_id)) %>%
  #   filter(!is.na(Terms) & !is.na(GDM_diag) & !is.na(Edu) & !is.na(Race) & !is.na(Nativity))
  filtered_demogra_data <- demogra_data %>%
    mutate(ppt_id = as.character(ppt_id))
  
  # Filter peak area data to include only participants with complete demographic data
  filtered_peak_area_long <- pos_peak_area_averaged %>%
    pivot_longer(-peak_id, names_to = "ppt_id", values_to = "peak_area") %>%
    left_join(filtered_demogra_data, by = "ppt_id")
  
  # Step 2: Define demographic variables for analysis
  # demographic_vars <- c("Terms", "GDM_diag", "Edu")
  # demographic_vars <- c("Terms", "GDM_diag", "Edu",'Race','Nativity')
  demographic_vars <- c("Terms")
  
  
  # Filepath for saving plots
  output_dir <- "D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/stats_test_positive/positive_112chemicals_20241229/"
  
  # Step 3: Hierarchical Clustering Function,, scale by participants
  # perform_hierarchical_clustering <- function(data, demo_var) {
  #   
  #   # Filter data for the specific demographic variable
  #   data_filtered <- data %>%
  #     filter(!is.na(.data[[demo_var]]))
  #   
  #   # Transform data into wide format for clustering
  #   data_wide <- data_filtered %>%
  #     select(peak_id, ppt_id, peak_area) %>%
  #     pivot_wider(names_from = ppt_id, values_from = peak_area) %>%
  #     column_to_rownames("peak_id")
  #   
  #   # Scale the data
  #   # scaled_data <- scale(data_wide, center = TRUE, scale = TRUE)
  #   scaled_data <- t(apply(data_wide, 1, function(row) scale(row, center = TRUE, scale = TRUE)))
  #   
  #   # Perform Hierarchical Clustering
  #   dist_matrix <- dist(scaled_data, method = "euclidean")  # Compute distance matrix
  #   hc <- hclust(dist_matrix, method = "complete")          # Perform hierarchical clustering
  #   
  #   # Cut the dendrogram into clusters
  #   clusters <- cutree(hc, k = 4)  # Cut into 4 clusters
  #   
  #   # Add cluster assignments to the data
  #   cluster_assignments <- data.frame(
  #     peak_id = rownames(scaled_data),
  #     cluster = as.factor(clusters)
  #   )
  #   
  #   # Prepare annotations using ppt_id
  #   demo_annotations <- data_filtered %>%
  #     distinct(ppt_id, Terms, GDM_diag, Edu,Race,Nativity) %>%
  #     column_to_rownames("ppt_id")
  #   
  #   # # Perform PCA
  #   # pca_result <- prcomp(scaled_data, scale. = TRUE)
  #   # pca_data <- as.data.frame(pca_result$x) %>%
  #   #   rownames_to_column("peak_id") %>%
  #   #   left_join(cluster_assignments, by = "peak_id")
  #   # 
  #   # # PCA Plot (without participant ID labels)
  #   # pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  #   #   geom_point(size = 4, alpha = 0.7) +
  #   #   labs(
  #   #     title = paste("PCA Clustering by", demo_var),
  #   #     x = "PC1 (Principal Component 1)", 
  #   #     y = "PC2 (Principal Component 2)",
  #   #     color = "Cluster"
  #   #   ) +
  #   #   scale_color_manual(values = c("red", "blue", "green", "purple")) +
  #   #   theme_minimal()
  #   # 
  #   # # Save the PCA plot
  #   # ggsave(
  #   #   filename = paste0(output_dir, "PCA_Plot_", demo_var, ".tif"),
  #   #   plot = pca_plot,
  #   #   width = 8,
  #   #   height = 6,
  #   #   dpi = 150
  #   # )
  #   
  #   print(pca_plot)
  #   
  #   # Heatmap Visualization
  #   heatmap_file <- paste0(output_dir, "Heatmap_", demo_var, ".png")
  #   png(heatmap_file, width = 1200, height = 1000, res = 150)
  #   pheatmap(
  #     scaled_data,
  #     cluster_rows = TRUE,            # Use hierarchical clustering for rows
  #     cluster_cols = TRUE,            # Use hierarchical clustering for columns
  #     annotation_col = demo_annotations, # Add demographic annotations
  #     main = paste("Hierarchical Clustering Heatmap -", demo_var),
  #     color = colorRampPalette(c("blue", "white", "red"))(50),
  #     fontsize_row = 6,
  #     show_colnames = FALSE  # Hide participant IDs at the bottom
  #   )
  #   
  #   dev.off()
  #   
  #   return(list(
  #     pca_plot = pca_plot,
  #     hc_result = hc,
  #     cluster_assignments = cluster_assignments
  #   ))
  # }
  
  
  # Step 4: Run hierarchical clustering for each demographic variable
  # results <- list()
  # 
  # for (var in demographic_vars) {
  #   message("Running hierarchical clustering for: ", var)
  #   results[[var]] <- perform_hierarchical_clustering(filtered_peak_area_long, var)
  # }


  
  ##scale by compound
  perform_hierarchical_clustering_bycmp <- function(data, demo_var) {
    
    # Filter data for the specific demographic variable
    data_filtered <- data %>%
      filter(!is.na(.data[[demo_var]]))
    
    # Transform data into wide format for clustering
    data_wide <- data_filtered %>%
      select(peak_id, ppt_id, peak_area) %>%
      pivot_wider(names_from = ppt_id, values_from = peak_area) %>%
      column_to_rownames("peak_id")
    
    # Scale the data by compound (row-wise)
    scaled_data <- t(apply(data_wide, 1, function(row) scale(row, center = TRUE, scale = TRUE)))
    
    # Perform Hierarchical Clustering
    dist_matrix <- dist(scaled_data, method = "euclidean")  # Compute distance matrix
    hc <- hclust(dist_matrix, method = "complete")          # Perform hierarchical clustering
    
    # Cut the dendrogram into clusters
    clusters <- cutree(hc, k = 2)  # Cut into 4 clusters
    
    # Add cluster assignments to the data
    cluster_assignments <- data.frame(
      peak_id = rownames(scaled_data),
      cluster = as.factor(clusters)
    )
    
    # Prepare annotations using ppt_id
    demo_annotations <- data_filtered %>%
      # distinct(ppt_id, Terms, GDM_diag, Edu, Race, Nativity) %>%
      distinct(ppt_id, Terms) %>%
      column_to_rownames("ppt_id")
    
    # Perform PCA
    # pca_result <- prcomp(scaled_data, scale. = TRUE)
    # pca_data <- as.data.frame(pca_result$x) %>%
    #   rownames_to_column("peak_id") %>%
    #   left_join(cluster_assignments, by = "peak_id")
    # 
    # # PCA Plot (without participant ID labels)
    # pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
    #   geom_point(size = 4, alpha = 0.7) +
    #   labs(
    #     title = paste("PCA Clustering by", demo_var),
    #     x = "PC1 (Principal Component 1)", 
    #     y = "PC2 (Principal Component 2)",
    #     color = "Cluster"
    #   ) +
    #   scale_color_manual(values = c("red", "blue", "green", "purple")) +
    #   theme_minimal()
    # 
    # # Save the PCA plot
    # ggsave(
    #   filename = paste0(output_dir, "PCA_Plot_", demo_var, ".tif"),
    #   plot = pca_plot,
    #   width = 8,
    #   height = 6,
    #   dpi = 150
    # )
    # 
    # print(pca_plot)
    
    # Heatmap Visualization
    heatmap_file <- paste0(output_dir, "Heatmap_", demo_var, ".png")
    png(heatmap_file, width = 1200, height = 1000, res = 150)
    pheatmap(
      scaled_data,
      cluster_rows = TRUE,            # Use hierarchical clustering for rows
      cluster_cols = TRUE,            # Use hierarchical clustering for columns
      annotation_col = demo_annotations, # Add demographic annotations
      main = paste("Hierarchical Clustering Heatmap -", demo_var),
      color = colorRampPalette(c("blue", "white", "red"))(50),
      fontsize_row = 6,
      show_colnames = FALSE  # Hide participant IDs at the bottom
    )
    
    dev.off()
    
    return(list(
      # pca_plot = pca_plot,
      hc_result = hc,
      cluster_assignments = cluster_assignments
    ))
  }
  
  results <- list()
  
  for (var in demographic_vars) {
    message("Running hierarchical clustering for: ", var)
    results[[var]] <- perform_hierarchical_clustering_bycmp(filtered_peak_area_long, var)
  }
  
  

  
  
  
