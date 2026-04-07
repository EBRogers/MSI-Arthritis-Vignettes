#' List all objects under a given prefix in a public S3 bucket.
#' Uses the S3 ListObjectsV2 REST API and parses the XML response with base R.
#' Handles pagination automatically for prefixes with >1000 objects.
s3_list_objects <- function(bucket_url, prefix) {
  all_keys  <- character(0)
  all_sizes <- numeric(0)
  continuation <- ""
  is_truncated <- TRUE
  
  while (is_truncated) {
    # Build the ListObjectsV2 request URL
    list_url <- paste0(bucket_url, "/?list-type=2&max-keys=1000&prefix=",
                       utils::URLencode(prefix, reserved = TRUE))
    if (nzchar(continuation)) {
      list_url <- paste0(list_url, "&continuation-token=",
                         utils::URLencode(continuation, reserved = TRUE))
    }
    
    # Fetch and read the XML response
    xml <- tryCatch(
      readLines(url(list_url), warn = FALSE),
      error = function(e) {
        warning("Failed to list objects for prefix '", prefix, "': ",
                conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(xml)) break
    xml <- paste(xml, collapse = "\n")
    
    # Extract <Key> values
    keys <- regmatches(xml, gregexpr("<Key>[^<]+</Key>", xml))[[1]]
    keys <- gsub("</?Key>", "", keys)
    
    # Extract <Size> values
    sizes <- regmatches(xml, gregexpr("<Size>[^<]+</Size>", xml))[[1]]
    sizes <- as.numeric(gsub("</?Size>", "", sizes))
    
    all_keys  <- c(all_keys, keys)
    all_sizes <- c(all_sizes, sizes)
    
    # Check for more pages
    trunc_match <- regmatches(xml,
                              regexpr("<IsTruncated>[^<]+</IsTruncated>", xml))
    is_truncated <- length(trunc_match) > 0 &&
      grepl("true", trunc_match, ignore.case = TRUE)
    
    if (is_truncated) {
      cont_match <- regmatches(
        xml, regexpr("<NextContinuationToken>[^<]+</NextContinuationToken>", xml)
      )
      continuation <- gsub("</?NextContinuationToken>", "", cont_match)
    }
  }
  
  # Filter out directory placeholders and .DS_Store files
  keep <- all_sizes > 0 &
    !grepl("/$", all_keys) &
    !grepl("\\.DS_Store$", all_keys)
  
  data.frame(key = all_keys[keep], size = all_sizes[keep],
             stringsAsFactors = FALSE)
}

#' Format bytes into a human-readable string.
human_size <- function(bytes) {
  if (bytes >= 1073741824) {
    sprintf("%.1f GB", bytes / 1073741824)
  } else if (bytes >= 1048576) {
    sprintf("%.1f MB", bytes / 1048576)
  } else if (bytes >= 1024) {
    sprintf("%.1f KB", bytes / 1024)
  } else {
    paste(bytes, "bytes")
  }
}

#' Download a single file from the public S3 bucket.
#' Skips files that already exist locally with the correct size.
#' Uses a .part temp file so interrupted downloads don't leave corrupt files.
s3_download_file <- function(bucket_url, key, dest_dir = ".") {
  local_path <- file.path(dest_dir, key)
  local_dir  <- dirname(local_path)
  remote_url <- paste0(bucket_url, "/", key)
  
  # Create the local directory tree if needed
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Skip if the file already exists and the size matches
  if (file.exists(local_path)) {
    local_size <- file.info(local_path)$size
    return(invisible(TRUE))
  }
  
  # Download to a temporary .part file first

  temp_path <- paste0(local_path, ".part")
  if (file.exists(temp_path)) file.remove(temp_path)
  
  success <- tryCatch({
    download.file(remote_url, destfile = temp_path, mode = "wb", quiet = TRUE)
    file.rename(temp_path, local_path)
    TRUE
  }, error = function(e) {
    warning("  Failed: ", key, " - ", conditionMessage(e))
    if (file.exists(temp_path)) file.remove(temp_path)
    FALSE
  })
  
  invisible(success)
}