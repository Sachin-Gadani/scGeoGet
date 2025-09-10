#' Download GEO Supplementary Files
#'
#' Downloads supplementary files from a GEO accession and organizes them
#'
#' @param accession Character. GEO accession number
#' @param target_dir Character. Directory to store downloaded files
#' @param verbose Logical. Print progress messages
#'
#' @return Character vector of downloaded file paths
#' @export
download_geo_data <- function(accession, target_dir, verbose = TRUE) {
  
  # Create target directory if it doesn't exist
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  # Create subdirectory for this accession
  accession_dir <- file.path(target_dir, accession)
  if (!dir.exists(accession_dir)) {
    dir.create(accession_dir, recursive = TRUE)
  }
  
  if (verbose) cat("Downloading files to:", accession_dir, "\n")
  
  tryCatch({
    # Download supplementary files using GEOquery
    file_info <- GEOquery::getGEOSuppFiles(
      GEO = accession,
      makeDirectory = FALSE,
      baseDir = accession_dir,
      fetch_files = TRUE
    )
    
    # Get the downloaded file paths
    downloaded_files <- rownames(file_info)
    
    if (length(downloaded_files) == 0) {
      stop("No supplementary files found for accession: ", accession)
    }
    
    # Check file integrity
    valid_files <- check_file_integrity(downloaded_files)
    if (!all(valid_files)) {
      warning("Some files appear corrupted or empty: ", 
              paste(downloaded_files[!valid_files], collapse = ", "))
      downloaded_files <- downloaded_files[valid_files]
    }
    
    if (verbose) {
      cat("Successfully downloaded", length(downloaded_files), "files:\n")
      for (f in downloaded_files) {
        cat("  -", basename(f), "\n")
      }
    }
    
    return(downloaded_files)
    
  }, error = function(e) {
    stop("Failed to download GEO data: ", e$message)
  })
}