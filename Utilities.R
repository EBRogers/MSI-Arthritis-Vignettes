#### Constants ####

SUBJECTS <- list() 
SUBJECTS[["CAD"]] <- c("6919", "6924", "6938", "6944") 
SUBJECTS[["OA"]] <- c("6891", "6908", "6949", "6954") 

CONDITION_IDS <- c("CAD", "OA") 
CONDITIONS <- c("Control", "Osteoarthritis") 

TISSUES <- c("lateral", "medial") 

name_tibble <- tribble(
  ~run,               ~subject, ~condition,        ~tissue,
  '6919-CAD-lateral',  1,       'Control',         'Lateral',
  '6919-CAD-medial',   1,       'Control',         'Medial',
  '6924-CAD-lateral',  2,       'Control',         'Lateral',
  '6924-CAD-medial',   2,       'Control',         'Medial',
  '6938-CAD-lateral',  3,       'Control',         'Lateral',
  '6938-CAD-medial',   3,       'Control',         'Medial',
  '6944-CAD-lateral',  4,       'Control',         'Lateral',
  '6944-CAD-medial',   4,       'Control',         'Medial',
  '6891-OA-lateral',   5,       'Osteoarthritis',  'Lateral',
  '6891-OA-medial',    5,       'Osteoarthritis',  'Medial',
  '6908-OA-lateral',   6,       'Osteoarthritis',  'Lateral',
  '6908-OA-medial',    6,       'Osteoarthritis',  'Medial',
  '6949-OA-lateral',   7,       'Osteoarthritis',  'Lateral',
  '6949-OA-medial',    7,       'Osteoarthritis',  'Medial',
  '6954-OA-lateral',   8,       'Osteoarthritis',  'Lateral',
  '6954-OA-medial',    8,       'Osteoarthritis',  'Medial'
) %>%
  mutate(tid = stringr::str_extract(run, "\\d{4}") %>% 
           paste0(ifelse(grepl("-CAD-", run), "-C", "-O"), 
                  ifelse(grepl("-lateral", run), "L", "M")),
         run = sub("^(([^-]*-){2}[^-]*)-.*$", "\\1", run),
         SUBJECT = stringr::str_extract(run, "\\d{4}"),
         CONDITION = ifelse(grepl("-CAD-", run), "CAD", "OA"),
         TISSUE = ifelse(grepl("-lateral", run), "lateral", "medial"))

#### Minor functions

flip_vertically <- function(m){
  max_y <- max(pData(m)$y)
  pData(m)$y <- max_y - pData(m)$y + 1
  return(m)
}

"%ni%" <- Negate("%in%")

is.na.durable <- function(x) {
  if (is.atomic(x)) {
    # Check if any element in the atomic vector is NA
    return(any(is.na(x)))
  } else if (is.list(x)) {
    return(FALSE)  # NA values within lists need to be handled individually
  } else {
    return(FALSE)  # For all other types of objects, return FALSE
  }
}

isolate_filename <- function(full_path) {
  parts <- strsplit(full_path, "/")[[1]]
  return(parts[[length(parts)]])
}

manage_directory <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, full.names = TRUE)
    file.remove(files)
    dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    unlink(dirs, recursive = TRUE)
  } else {
    dir.create(path, recursive = TRUE)
  }
}

#### MSE write functions ####

generate_dir_structure <- function(root) {
  # create root directory if it doesn't exist 
  if (!dir.exists(root)) {
    dir.create(root, recursive=TRUE)
  } 
  
  # CAD and OA 
  for (condition in CONDITION_IDS) {
    cond_dir <- file.path(root, condition) 
    dir.create(cond_dir, showWarnings=FALSE) 
    
    # each subject in CAD and OA 
    for (subject in SUBJECTS[[condition]]) {
      subject_dir <- file.path(cond_dir, subject) 
      dir.create(subject_dir, showWarnings=FALSE) 
      
      # for each tissue 
      for (tissue in TISSUES) {
        tissue_dir <- file.path(subject_dir, tissue) 
        dir.create(tissue_dir, showWarnings=FALSE)
      }
    }
  } 
} 

write_MSI_list <- function(mse_list, root) { 
  generate_dir_structure(root = root)
  for (i in c(1:length(mse_list))) {
    name_tibble_row <- name_tibble[name_tibble$run == runNames(mse_list[[i]]),]
    subject <- name_tibble_row$SUBJECT
    condition <- name_tibble_row$CONDITION
    tissue <- name_tibble_row$TISSUE
    
    PREFIX <- paste0(subject, "-", condition, 
                     "-", tissue, "-") 
    
    writeMSIData(
      mse_list[[i]], 
      file.path(root, condition, subject, tissue, 
                paste0(name_tibble_row$run, ".imzML")) 
    )
  }
}

#### Clustering functions ####

aggregate_feature_groups <- function(matrix_data, group_vector, crawl_along, agg_fun = "mean", fill.na=NULL,na.rm=FALSE) {
  
  stopifnot(agg_fun %in% c("mean", "median"))
  
  if (!is.null(fill.na)) {
    matrix_data <- replace(matrix_data, is.na(matrix_data), fill.na)
  }
  # Create a list of row indices for each group
  unique_groups <- crawl_along
  # Pre-allocate a result matrix
  result <- matrix(nrow = length(unique_groups),
                   ncol = ncol(matrix_data))
  for (i in seq_along(unique_groups)) {
    # Find the rows belonging to the current group
    group_rows <- which(group_vector == unique_groups[i])
    # Compute the column-wise mean of the selected rows
    if (agg_fun == "mean") {
      result[i, ] <- base::colMeans(matrix_data[group_rows, , drop = FALSE], na.rm = na.rm)
    } else if (agg_fun == "median") {
      result[i, ] <- matrixStats::colMedians(matrix_data[group_rows, , drop = FALSE], na.rm = na.rm) 
    }
  }
  return(result)
}

condense_groups <- function(mse, group_frame, agg_fun = "mean", fill.na=NULL, na.rm=FALSE, verbose = FALSE) {
  
  stopifnot(agg_fun %in% c("mean", "median"))
  
  stopifnot("group_id" %in% names(group_frame) & 
              "mz" %in% names(group_frame) & 
              "feature_id" %in% names(group_frame))
  
  fd <- merge(
    tibble(
      mz = featureData(mse)$mz,
      feature_id = 1:length(featureData(mse)$mz)
    ),
    group_frame,
    by = "feature_id",
    suffixes = c("", "_group_frame")
  )
  
  tol <- 1e-10
  
  stopifnot(
    fd %>%
      mutate(check = mz == mz_group_frame,
             delta = abs(mz - mz_group_frame) < tol) %>%
      group_by(delta) %>%
      summarise(n = n()) %>%
      pull(n) == c(nrow(group_frame))
  )
  
  if (verbose) {
    print("Beginning matter_mat grouping")
  }
  
  if (class(intensity(mse))[1] == "mat_matrix") {
    dn <- intensity(mse)@dimnames
  } else {
    dn <- NULL
  }
  
  condensed_fd <- fd %>%
    group_by(group_id) %>%
    summarise(mz = mean(mz),
              group_size = n(),
              group = group_size > 1) %>%
    ungroup() %>%
    arrange(mz) %>%
    mutate(new_feature_id = row_number())
  
  new_intensity_matter_mat <- matter::matter_mat(aggregate_feature_groups(intensity(mse), 
                                                                          fd$group_id,
                                                                          condensed_fd$group_id,
                                                                          agg_fun = agg_fun,
                                                                          fill.na=fill.na,na.rm=na.rm), 
                                                 dimnames = dn)
  
  mdf <- MassDataFrame(
    mz = condensed_fd$mz,
    group_size = condensed_fd$group_size,
    group = condensed_fd$group
  )
  
  intensity(mse) <- new_intensity_matter_mat
  
  featureData(mse) <- mdf
  
  return(mse)
}

#### Modeling Functions ####
get_convergence_code = function(m) {
  res <- tryCatch({
    m@optinfo$conv$lme4$code
  }, error = function(e) {
    NA
  })
  if (is.null(res)) {
    res <- NA
  }
  return(res)
}
