
SUBJECTS <- list() 
SUBJECTS[["CAD"]] <- c("6919", "6924", "6938", "6944") 
SUBJECTS[["OA"]] <- c("6891", "6908", "6949", "6954") 

CONDITION_IDS <- c("CAD", "OA") 
CONDITIONS <- c("Control", "Osteoarthritis") 

SUBTISSUES <- c("lateral", "medial") 
PROC_MILESTONE_IMPORT <- "recalibrated" 

generate_run_names <- function(
    subjects, condition_code, proc_milestone) { 
    run <- sapply(subjects, function(subject) {
        prefix <- paste0(subject, "-", condition_code, "-") 
        
        out <- c()
        for (subtissue in SUBTISSUES) {
            out <- c(out, 
                paste0(prefix, subtissue, "-", proc_milestone)) 
        }
        out
    }) 
    run |> as.vector() 
} 

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
            
            # for each subtissue 
            for (subtissue in SUBTISSUES) {
                subtissue_dir <- file.path(subject_dir, subtissue) 
                dir.create(subtissue_dir, showWarnings=FALSE)
            }
        }
    } 
} 

rename_runnames <- function(mse_list, suffix) {
    for (i in c(1:length(mse_list))) {
        runname <- runNames(mse_list[[i]]) 
        runNames(mse_list[[i]]) <- paste(
            c(
                strsplit(runname, split="-")[[1]][1:3], 
                suffix
            ), 
            collapse = "-"
        )
    } 
    mse_list
} 

write_MSI_list <- function(mse_list, root) { 
    generate_dir_structure(root = root)
    names <- sapply(mse_list, runNames) |> unname() 
    for (i in c(1:length(names))) {
        subject <- EXP_META$subject[i] 
        condition <- EXP_META$condition[i] 
        subtissue <- EXP_META$subtissue[i] 
        
        PREFIX <- paste0(subject, "-", condition, 
            "-", subtissue, "-") 
        
        writeMSIData(
            mse_list[[i]], 
            file.path(root, condition, subject, subtissue, 
                      paste0(names[i], ".imzML")) 
        )
    }
}


EXP_META <- tibble(
    run = c(
        generate_run_names(SUBJECTS$CAD, "CAD", 
            PROC_MILESTONE_IMPORT), 
        generate_run_names(SUBJECTS$OA, "OA", 
            PROC_MILESTONE_IMPORT)
    ), 
    subject = c(
        rep(SUBJECTS$CAD, 2) |> matrix(ncol=2) |> t() |> as.vector(), 
        rep(SUBJECTS$OA, 2) |> matrix(ncol=2) |> t() |> as.vector() 
    ), 
    condition = c(rep("CAD", 8), rep("OA", 8)), 
    subtissue = c(rep(SUBTISSUES, 8))
) 

EXP_META <- EXP_META |> 
      mutate(tid = stringr::str_extract(run, "\\d{4}") |> 
           paste0(ifelse(grepl("-CAD-", run), "-C", "-O"), 
                  ifelse(grepl("-lateral-", run), "L", "M"))) 


MSE_FILENAMES <- list.files(
    file.path(DATAPATH, CURRPATH), 
    pattern = "\\.imzML$", 
    recursive = TRUE
) 

isolate_filename <- function(full_path) {
  return(basename(full_path))
}

flip_vertically <- function(m){
  max_y <- max(pData(m)$y)
  pData(m)$y <- max_y - pData(m)$y
  return(m)
} 

# MSE_FILENAMES <- MSE_FILENAMES[-1] 
mses = list()
for (f in MSE_FILENAMES) {
  just_fname <- isolate_filename(f)
  m <- readMSIData(file.path(DATAPATH, CURRPATH, f),verbose=F) 
  mses[[just_fname]] <- m 
} 


medians_by_pixel <- function(mse) {
    mse |> 
        spectra() |> 
        as.matrix(nrow=dim(mse)[1]) |> 
        matrixStats::colMedians(na.rm=TRUE) 
} 

convert_zeros_to_NA <- function(vector_ints) {
    vector_ints[vector_ints == 0.0] <- NA 
    vector_ints
} 

op <- function(vector_ints) { 
    thresh <- quantile(vector_ints, 0.95) 
    if (thresh == 0.0) thresh <- max(vector_ints) 
    
    vector_ints |> 
        pmin(thresh) |> 
        convert_zeros_to_NA() 
}

clip_intensities <- function(mse) { 
    .mse <- mse |> spectra() |> as.matrix()
    
    clipped_spectra <- apply(
            .mse, 
            MARGIN = 2, 
            op 
        )
    
    intensity(mse) <- clipped_spectra
    mse
} 


NAize_sparse_spectra <- function(mse) { 
    spa <- mse |> spectra() 
    nf <- dim(mse)[1] 
    np <- dim(mse)[2]

    for (i in c(1:np)) {
        ct <- (spa[, i] > 2) |> which() |> length() 
        if (ct < nf/2) spa[, i] <- NA 
    }

    spectra(mse) <- spa 
    mse
} 

normalize_with_median <- function(mse) {
    scaled_spectra <- mse |> 
        spectra() |> 
        apply(
            MARGIN = 2, 
            function(x) {
                med <- median(x, na.rm=TRUE) 
                
                scaler <- MEDIAN / med 
                x <- scaler * x 
            }
        ) 
    intensity(mse) <- scaled_spectra 
    mse 
} 


.identical <- function(vec, target) {
    vec <- vec[!is.na(vec)] 
    all.equal(vec, rep(target, length(vec)))
}





















