
scale_row_sd_mean_subset <- function(matrix, 
                                     target_sd = 90, 
                                     target_mean = NULL, 
                                     columns_to_scale = rep(TRUE, ncol(matrix)), 
                                     verbose = FALSE) {
  if (length(columns_to_scale) != ncol(matrix)) {
    stop("columns_to_scale must be a logical vector with length equal to the number of columns in the matrix")
  }
  
  if (length(target_sd) == 1) {
    target_sds <- rep(target_sd, nrow(matrix))
  } else if (length(target_sd) == nrow(matrix)) {
    target_sds <- target_sd
  } else {
    stop("Length of target_sd must be either 1 or equal to the number of rows in the matrix")
  }
  
  scaled_matrix <- matrix
  
  row_sds <- apply(matrix[, columns_to_scale, drop = FALSE], 1, sd, na.rm = TRUE)
  row_means <- apply(matrix[, columns_to_scale, drop = FALSE], 1, mean, na.rm = TRUE)
  scaling_factors <- target_sds / row_sds
  
  if (is.null(target_mean)) {
    target_means <- row_means
  } else if (length(target_mean) == 1) {
    if (target_mean <= 0) stop("Target mean must be positive")
    target_means <- rep(target_mean, nrow(matrix))
  } else if (length(target_mean) == nrow(matrix)) {
    if (any(target_mean <= 0)) stop("All target means must be positive")
    target_means <- target_mean
  } else {
    stop("Length of target_mean must be either 1 or equal to the number of rows in the matrix")
  }
  
  for (i in 1:nrow(matrix)) {
    current_row <- matrix[i, columns_to_scale]
    na_positions <- is.na(current_row)
    
    if (all(na_positions)) next
    
    non_na_values <- current_row[!na_positions]
    centered_non_na <- non_na_values - row_means[i]
    scaled_non_na <- centered_non_na * scaling_factors[i]
    adjusted_non_na <- scaled_non_na + target_means[i]
    
    if (any(adjusted_non_na < 0)) {
      min_centered_original <- min(centered_non_na)
      if (min_centered_original < 0) {
        reduced_scaling_factor <- target_means[i] / abs(min_centered_original)
        actual_scaling_factor <- min(scaling_factors[i], reduced_scaling_factor)
        scaled_non_na <- centered_non_na * actual_scaling_factor
        adjusted_non_na <- scaled_non_na + target_means[i]
        actual_sd <- sd(adjusted_non_na)
        if (verbose) {
          cat("Row", i, "SD reduced from", target_sds[i], "to", actual_sd, "to avoid negative values\n")
        }
      }
    }
    
    scaled_matrix[i, which(columns_to_scale)[!na_positions]] <- adjusted_non_na
  }
  
  return(scaled_matrix)
}




run_models_no_me <- function(mse_frame, formu = "mean ~ condition", verbose = FALSE) {
  out <- list()
  model_results <-NULL
  models <- list()
  for (f in unique(mse_frame$feature_id)) {
    if (verbose) print(paste0("Feature ",f))
    if (inherits(mse_frame, "data.table")) {
      # Filter using data.table syntax
      if (verbose) print("using dt filtering")
      data <- mse_frame[feature_id == f]
    } else {
      # Filter using dplyr syntax for tibble/data.frame
      data <- mse_frame %>% filter(feature_id == f)
    }
    m <- NULL
    m <- lm(formu, data = data)
    interm <- rbind(tidy(m)
                    %>% mutate(feature_id = f,
                               model_id = "means"
                               #total_var = sum(as.data.frame(VarCorr(m))$vcov),
                    )
    )
    
    model_results <- rbind(model_results, 
                           interm %>% pivot_wider(id_cols = c(feature_id, model_id), 
                                                  names_from = term, 
                                                  values_from = c(estimate, std.error, statistic, p.value)))
    
    models[[f]] <- m
    if (inherits(mse_frame, "data.table")) {
      # Filter using data.table syntax
      data <- mse_frame[feature_id == f]
    }
  }
  
  model_results_long <- model_results %>%
    pivot_longer(
      cols = starts_with("estimate_") | starts_with("std.error_") | starts_with("statistic_") | starts_with("df_") | starts_with("p.value_"),
      names_to = c("model_item", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    )
  
  out[["models"]] <- models
  out[["model_results"]] <- model_results
  out[["model_results_long"]] <- model_results_long
  
  return(out)
  
}

autoscore_sscl <- function(ssc_obj, mse_obj){
  class_tib <- tibble(x = pData(mse_obj)$x, 
                      y = pData(mse_obj)$y, 
                      oracle_roi = pData(mse_obj)$circle_2)
  rf <- NULL
  for (id in names(ssc_obj)) {
    rf <- bind_rows(rf,
                    class_tib %>% 
                      mutate(class = ssc_obj[[id]]$class) %>% 
                      summarise(n = n(), .by = c(class, oracle_roi)) %>% 
                      pivot_wider(names_from = oracle_roi, id_cols = class, values_from = n) %>% 
                      mutate(across(c(`FALSE`,`TRUE`), ~replace_na(.x, 0)),
                             score = (`FALSE`/-6232) + (`TRUE`/3768)) %>% 
                      select(class, score) %>%
                      mutate(id = id))
  }
  return(rf)
}

autoscore_skml <- function(skm_obj, mse_obj){
  class_tib <- tibble(x = pData(mse_obj)$x, 
                      y = pData(mse_obj)$y, 
                      oracle_roi = pData(mse_obj)$circle_2)
  
  class_tib %>% 
    mutate(class = skm_obj$cluster) %>% 
    summarise(n = n(), .by = c(class, oracle_roi)) %>% 
    pivot_wider(names_from = oracle_roi, id_cols = class, values_from = n) %>% 
    mutate(across(c(`FALSE`,`TRUE`), ~replace_na(.x, 0)),
           score = (`FALSE`/-6232) + (`TRUE`/3768)) %>% 
    select(class, score)
}

analyze_spatial_clusters <- function(sscl, 
                                     sim_list_adj, 
                                     samples, 
                                     n_top_per_condition = 4,
                                     n_top_features = 20,
                                     output_dir = FALSE) {
  #' Analyze spatial shrunken centroids clustering results
  #'
  #' @param sscl List of spatialShrunkenCentroids objects from Cardinal package
  #' @param sim_list_adj List of MSImagingExperiment objects from Cardinal package  
  #' @param samples Data frame with columns: subject, subtissue, condition, id
  #' @param n_top_per_condition Number of top scoring clusters to select per condition (default: 4)
  #' @param n_top_features Number of top features to extract per cluster (default: 20)
  #' @param output_dir Directory name for saving images (default: "ssim_temp")
  #'
  #' @return List containing: selected_clusters, sim_to_test, samples_to_test, unique_top_20_all_samples
  
  # Score and filter clusters
  score_filtered_sscls <- bind_rows(lapply(1:16, function(i){
    autoscore_sscl(sscl[[i]], sim_list_adj[[i]])
  }), .id = "sim_id") %>% 
    filter(score == max(score), .by = sim_id) %>% 
    arrange(as.numeric(sim_id)) %>% 
    bind_cols(samples %>% select(!id)) %>% 
    arrange(desc(score)) %>% 
    slice_head(n = n_top_per_condition, by = condition)
  
  # Create selected clusters list
  selected_clusters <- list()
  for (r in 1:nrow(score_filtered_sscls)) {
    selected_clusters[[r]] <- list("sim_adj_id" = as.integer(score_filtered_sscls[[r,"sim_id"]]),  
                                   "ssc_index" = score_filtered_sscls[[r,"id"]], 
                                   "roi_cluster_id" = score_filtered_sscls[[r,"class"]])
  }
  
  # Set up output directory
  if (output_dir != FALSE) {
  manage_directory(output_dir)
  }
  # Process selected clusters
  sim_to_test <- list()
  samples_to_test <- NULL
  unique_top_20_all_samples <- c()
  i = 1
  
  for (sc in selected_clusters) {
    # Save cluster image
    if (output_dir != FALSE) {
    png(filename = paste0(output_dir, "/", sc[[1]], "_", sc[[2]], "_", sc[[3]], "_ssc.png"))
    print(image(sscl[[sc$sim_adj_id]][[sc$ssc_index]], "class"))
    dev.off()
    }
    
    # Add ROI classification to pData
    pData(sim_list_adj[[sc$sim_adj_id]])$roi_ssc <- ifelse(sscl[[sc$sim_adj_id]][[sc$ssc_index]]$class == as.character(sc$roi_cluster_id),
                                                           "ssc_roi",
                                                           "other")
    m <- sim_list_adj[[sc$sim_adj_id]]
    
    # Summarize features
    sim_to_test[[i]] <- summarizeFeatures(m, 
                                          stat = c("mean", "sd"), 
                                          na.rm = T, 
                                          verbose = F, 
                                          groups = m$roi_ssc)
    
    # Collect top features
    unique_top_20_all_samples <- unique(c(unique_top_20_all_samples, 
                                          topFeatures(sscl[[sc$sim_adj_id]][[sc$ssc_index]], 
                                                      n = n_top_features)$i))
    
    # Collect samples
    samples_to_test <- bind_rows(samples_to_test, samples[sc$sim_adj_id,])
    i <- i + 1
    
    # Save image classification
    if (output_dir != FALSE) {
    png(filename = paste0(output_dir, "/", sc[[1]], "_", sc[[2]], "_", sc[[3]], "_imageclass.png"))
    print(image(sim_list_adj[[sc$sim_adj_id]], "roi_ssc"))
    dev.off()
    }
  }
  
  # Return all results as a list
  return(list(
    selected_clusters = selected_clusters,
    sim_to_test = sim_to_test,
    samples_to_test = samples_to_test,
    unique_top_20_all_samples = unique_top_20_all_samples
  ))
}

analyze_spatial_clusters_skm <- function(skml, 
                                     sim_list_adj, 
                                     samples, 
                                     n_top_per_condition = 4,
                                     n_top_features = 20,
                                     output_dir = FALSE) {
  #' Analyze spatial shrunken centroids clustering results
  #'
  #' @param skml List of spatialKMeans objects from Cardinal package
  #' @param sim_list_adj List of MSImagingExperiment objects from Cardinal package  
  #' @param samples Data frame with columns: subject, subtissue, condition, id
  #' @param n_top_per_condition Number of top scoring clusters to select per condition (default: 4)
  #' @param n_top_features Number of top features to extract per cluster (default: 20)
  #' @param output_dir Directory name for saving images (default: "ssim_temp")
  #'
  #' @return List containing: selected_clusters, sim_to_test, samples_to_test, unique_top_20_all_samples
  
  # Score and filter clusters
  score_filtered_skmls <- bind_rows(lapply(1:16, function(i){
    autoscore_skml(skml[[i]], sim_list_adj[[i]])
  }), .id = "sim_id") %>% 
    filter(score == max(score), .by = sim_id) %>% 
    arrange(as.numeric(sim_id)) %>% 
    bind_cols(samples %>% select(!id)) %>% 
    arrange(desc(score)) %>% 
    slice_head(n = n_top_per_condition, by = condition)
  
  # Create selected clusters list
  selected_clusters <- list()
  for (r in 1:nrow(score_filtered_skmls)) {
    selected_clusters[[r]] <- list("sim_adj_id" = as.integer(score_filtered_skmls[[r,"sim_id"]]), 
                                   "roi_cluster_id" = score_filtered_skmls[[r,"class"]])
  }
  
  # Set up output directory
  if (output_dir != FALSE) {
  manage_directory(output_dir)
  }
  
  # Process selected clusters
  sim_to_test <- list()
  samples_to_test <- NULL
  unique_top_20_all_samples <- c()
  i = 1
  
  for (sc in selected_clusters) {
    # Save cluster image
    if (output_dir != FALSE) {
      png(filename = paste0(output_dir, "/", sc[[1]], "_", sc[[2]], "_ssc.png"))
      print(image(skml[[sc$sim_adj_id]], "cluster"))
      dev.off()
    }
   
    
    # Add ROI classification to pData
    pData(sim_list_adj[[sc$sim_adj_id]])$roi_skm <- ifelse(as.integer(skml[[sc$sim_adj_id]]$cluster) == as.integer(sc$roi_cluster_id),
                                                           "skm_roi",
                                                           "other")
    m <- sim_list_adj[[sc$sim_adj_id]]
    
    # Summarize features
    sim_to_test[[i]] <- summarizeFeatures(m, 
                                          stat = c("mean", "sd"), 
                                          na.rm = T, 
                                          verbose = F, 
                                          groups = m$roi_skm)
    
    # Collect top features
    unique_top_20_all_samples <- unique(c(unique_top_20_all_samples, 
                                          topFeatures(skml[[sc$sim_adj_id]], 
                                                      n = n_top_features)$i))
    
    # Collect samples
    samples_to_test <- bind_rows(samples_to_test, samples[sc$sim_adj_id,])
    i <- i + 1
    
    # Save image classification
    if (output_dir != FALSE) {
      png(filename = paste0(output_dir, "/", sc[[1]], "_", sc[[2]], "_imageclass.png"))
      print(image(sim_list_adj[[sc$sim_adj_id]], "roi_skm"))
      dev.off()
    }
  }
  
  # Return all results as a list
  return(list(
    selected_clusters = selected_clusters,
    sim_to_test = sim_to_test,
    samples_to_test = samples_to_test,
    unique_top_20_all_samples = unique_top_20_all_samples
  ))
}



#### FOR SIM 2


scale_row_sd_mean <- function(matrix, target_sd = 90, target_mean = NULL, verbose = F) {
  if (length(target_sd) == 1) {
    target_sds <- rep(target_sd, nrow(matrix))
  } else if (length(target_sd) == nrow(matrix)) {
    target_sds <- target_sd
  } else {
    stop("Length of target_sd must be either 1 or equal to the number of rows in the matrix")
  }
  
  row_sds <- apply(matrix, 1, sd, na.rm = TRUE)
  row_means <- apply(matrix, 1, mean, na.rm = TRUE)
  scaling_factors <- target_sds / row_sds
  
  if (is.null(target_mean)) {
    target_means <- row_means
  } else if (length(target_mean) == 1) {
    if (target_mean <= 0) stop("Target mean must be positive")
    target_means <- rep(target_mean, nrow(matrix))
  } else if (length(target_mean) == nrow(matrix)) {
    if (any(target_mean <= 0)) stop("All target means must be positive")
    target_means <- target_mean
  } else {
    stop("Length of target_mean must be either 1 or equal to the number of rows in the matrix")
  }
  
  scaled_matrix <- matrix
  for (i in 1:nrow(matrix)) {
    current_row <- matrix[i, ]
    na_positions <- is.na(current_row)
    
    if (all(na_positions)) next
    
    non_na_values <- current_row[!na_positions]
    centered_non_na <- non_na_values - row_means[i]
    scaled_non_na <- centered_non_na * scaling_factors[i]
    adjusted_non_na <- scaled_non_na + target_means[i]
    
    if (any(adjusted_non_na < 0)) {
      min_centered_original <- min(centered_non_na)
      
      if (min_centered_original < 0) {
        reduced_scaling_factor <- target_means[i] / abs(min_centered_original)
        actual_scaling_factor <- min(scaling_factors[i], reduced_scaling_factor)
        scaled_non_na <- centered_non_na * actual_scaling_factor
        adjusted_non_na <- scaled_non_na + target_means[i]
        actual_sd <- sd(adjusted_non_na)
        if (verbose) {
          cat("Row", i, "SD reduced from", target_sds[i], "to", actual_sd, "to avoid negative values\n")
        }
      }
    }
    
    scaled_matrix[i, !na_positions] <- adjusted_non_na
  }
  
  return(scaled_matrix)
}



run_models <- function(mse_frame, formu = "mean ~ condition * subtissue + (1|subject)", verbose = FALSE ) {
  out <- list()
  model_results <-NULL
  models <- list()
  for (f in unique(mse_frame$feature_id)) {
    if (verbose) print(paste0("Feature ",f))
    if (inherits(mse_frame, "data.table")) {
      # Filter using data.table syntax
      if (verbose) print("using dt filtering")
      data <- mse_frame[feature_id == f]
    } else {
      # Filter using dplyr syntax for tibble/data.frame
      data <- mse_frame %>% filter(feature_id == f)
    }
    m <- NULL
    m <- lmer(formu, data = data)
    interm <- rbind(tidy(m, "fixed")
                    %>% mutate(feature_id = f,
                               model_id = "means",
                               #icc_adj = get_icc(m, "ICC_adjusted"),
                               #icc_unadj = get_icc(m, "ICC_unadjusted"),
                               icc = as.data.frame(VarCorr(m))$vcov[1]/sum(as.data.frame(VarCorr(m))$vcov),
                               total_var = sum(as.data.frame(VarCorr(m))$vcov),
                               #r2.marginal = get_gof(m)$r2.marginal,
                               isSingular = isSingular(m)
                               #convergenceCode = get_convergence_code(m))
                               )
                    )
    
    model_results <- rbind(model_results, 
                           interm %>% pivot_wider(id_cols = c(feature_id, model_id, icc, total_var, isSingular), 
                                                  names_from = term, 
                                                  values_from = c(estimate, std.error, statistic, df, p.value)))
    
    models[[f]] <- m
    if (inherits(mse_frame, "data.table")) {
      # Filter using data.table syntax
      data <- mse_frame[feature_id == f]
    }
  }
  
  model_results_long <- model_results %>%
    pivot_longer(
      cols = starts_with("estimate_") | starts_with("std.error_") | starts_with("statistic_") | starts_with("df_") | starts_with("p.value_"),
      names_to = c("model_item", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    )
  
  out[["models"]] <- models
  out[["model_results"]] <- model_results
  out[["model_results_long"]] <- model_results_long
  
  return(out)
  
}

run_models_MEMORY_EFF <- function(mse_frame, formu = "mean ~ condition * subtissue + (1|subject)", batch_size = 50, verbose = FALSE) {
  out <- list()
  model_results <- NULL
  models <- list()
  
  # Get unique feature IDs and determine number of batches
  features <- unique(mse_frame$feature_id)
  total_features <- length(features)
  num_batches <- ceiling(total_features / batch_size)
  
  for (batch in 1:num_batches) {
    # Determine start and end indices for current batch
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, total_features)
    batch_features <- features[start_idx:end_idx]
    
    if (verbose) {
      cat(sprintf("Processing batch %d of %d (features indicies %d to %d)\n", 
                  batch, num_batches, start_idx, end_idx))
    }
    
    for (f in batch_features) {
      if (verbose) cat(sprintf("Feature %s\n", f))
      
      # Filter data - using data.table approach when possible for efficiency
      if (inherits(mse_frame, "data.table")) {
        data <- mse_frame[feature_id == f]
      } else {
        data <- mse_frame %>% filter(feature_id == f)
      }
      
      # Run model
      tryCatch({
        m <- lmer(formu, data = data)
        
        # Extract model results
        interm <- tidy(m, "fixed") %>% 
          mutate(
            feature_id = f,
            model_id = "means",
            icc = as.data.frame(VarCorr(m))$vcov[1]/sum(as.data.frame(VarCorr(m))$vcov),
            total_var = sum(as.data.frame(VarCorr(m))$vcov),
            isSingular = isSingular(m)
          )
        
        # Add to results
        model_results <- rbind(model_results, 
                               interm %>% pivot_wider(
                                 id_cols = c(feature_id, model_id, icc, total_var, isSingular), 
                                 names_from = term, 
                                 values_from = c(estimate, std.error, statistic, df, p.value)
                               ))
        
        # Store only essential model information if needed
        models[[f]] <- m
        
      }, error = function(e) {
        if (verbose) {
          cat(sprintf("Error in feature %s: %s\n", f, e$message))
        }
      })
      
      # Clear some memory
      data <- NULL
      gc()
    }
    
    # Perform thorough garbage collection between batches
    gc(full = TRUE, verbose = verbose)
    
    # Save intermediate results to disk if desired
    # saveRDS(model_results, paste0("model_results_batch_", batch, ".rds"))
  }
  
  # Process final results
  model_results_long <- model_results %>%
    pivot_longer(
      cols = starts_with("estimate_") | starts_with("std.error_") | 
        starts_with("statistic_") | starts_with("df_") | starts_with("p.value_"),
      names_to = c("model_item", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    )
  
  out[["models"]] <- models
  out[["model_results"]] <- model_results
  out[["model_results_long"]] <- model_results_long
  
  return(out)
}

get_comparison_plotting_frame <- function(mse_frame, model_res, nmax = 150) {
  ttest_frame_condition <- mse_frame %>% 
    arrange(feature_id) %>%
    filter(subtissue == "1") %>%
    group_by(feature_id) %>% 
    group_modify(~tidy(t.test(mean ~ condition, data = .x, var.equal = T))) %>%
    ungroup() %>%
    mutate(SE = estimate / statistic,
           adj_p.value = p.adjust(p.value, method = "fdr"),
           old_feature_id = feature_id,
           feature_id = row_number())
  
  
  ttest_frame_subtissue <- mse_frame %>% 
    filter(condition == "1") %>%
    group_by(feature_id) %>% 
    arrange(subtissue, subject) %>%
    group_modify(~tidy(t.test(.x$mean[.x$subtissue == "0"],
                              .x$mean[.x$subtissue == "1"],
                              paired = T))) %>%
    ungroup() %>%
    mutate(SE = estimate / statistic,
           adj_p.value = p.adjust(p.value, method = "fdr"),
           old_feature_id = feature_id,
           feature_id = row_number())
  
  contrasts_means <- bind_rows(lapply(model_res, FUN = function(m){
    bind_rows(tidy(contrast(emmeans(m, ~ condition | subtissue, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="subtissue"),
              tidy(contrast(emmeans(m, ~ subtissue | condition, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="condition")) 
  })) %>% 
    group_by(term, subgroup) %>%
    mutate(feature_id = row_number(),
           adj_p.value = p.adjust(p.value, method = "fdr")) %>% ungroup()
  
  tests <- bind_rows(
    
    contrasts_means%>%
      filter(subgroup=="1") %>%
      mutate(significant=`p.value`< 0.05,
             truth = ifelse(as.numeric(feature_id) > nmax,
                            TRUE,
                            FALSE),
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding, term, .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>% 
      mutate(model = "Mixed Effects Model with All Variables",
             test = "T Test of Conditional Means",
             null = case_match(term,
                               "condition"~"In Medial:\nControl vs Osteoarthritis",
                               "subtissue"~"In Osteoarthritis:\nLateral vs Medial"),
             term = case_match(term,
                               "condition"~"Condition",
                               "condition:subtissue"~"Interaction",
                               "subtissue"~"Subsection")) %>%
      dplyr::select(finding, n, model, test, null, term),
    
    ttest_frame_condition %>% 
      mutate(significant = `p.value`< 0.05,
             truth = ifelse(as.numeric(feature_id) > nmax,
                            TRUE,
                            FALSE),
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding, .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>% 
      mutate(model = "Linear Model with Condition",
             null = "In Medial:\nControl vs Osteoarthritis",
             test = "Two Sample T Test",
             term = "Condition") %>%
      dplyr::select(finding, n, model, test, null, term),
    
    ttest_frame_subtissue %>% 
      mutate(significant=`p.value`< 0.05,
             truth = ifelse(as.numeric(feature_id) > nmax,
                            TRUE,
                            FALSE),
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding,  .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>%
      mutate(model = "Linear Model of Paired Differences\nBetween Subtissues",
             null = "In Osteoarthritis:\nLateral vs Medial",
             test = "Paired T Test",
             term = "Subsection") %>%
      dplyr::select(finding, n, model, test, null, term)
  ) %>% mutate(termc = case_match(term,
                                  "Condition"~"Control vs OA",
                                  "Interaction"~"Interaction",
                                  "Subsection"~"Lateral vs Medial")) %>%
    mutate(termc = ifelse(test == "F Test of Mean Squares",
                          term, 
                          termc)) %>%
    mutate(termc = factor(termc, 
                          levels = c("Condition",
                                     "Subsection",
                                     "Interaction",
                                     "Control vs OA",
                                     "Lateral vs Medial"),
                          labels = c("Condition\nAny Difference",
                                     "Subsection\nAny Difference",
                                     "Interaction\nAny Present",
                                     "Condition\nControl vs OA",
                                     "Subsection\nLateral vs Medial")))
  return(tests)
  
}

#SHARED

extract_data_for_modeling_alt <- function(sim_list, samples_df) {
  mse_frame <- NULL
  for (i in 1:length(sim_list)) {
    m <- sim_list[[i]]
    t <- as_tibble(fData(m)) %>% mutate(feature_id = 1:length(mz(m)),
                                        subject = samples_df[i,]$subject,
                                        condition = samples_df[i,]$condition,
                                        subtissue = samples_df[i,]$subtissue
    )
    mse_frame <- bind_rows(mse_frame,t)
  }
  mse_frame <- mse_frame %>% mutate(across(c(subject, condition, subtissue), as.factor))
  return(mse_frame)
}