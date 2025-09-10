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
    
    # Extract tar files if present
    extracted_files <- c()
    for (file in downloaded_files) {
      if (grepl("\\.tar$", file, ignore.case = TRUE)) {
        if (verbose) cat("Extracting tar archive:", basename(file), "\n")
        
        # Extract tar file in the same directory
        utils::untar(file, exdir = dirname(file))
        
        # List extracted files
        tar_contents <- utils::untar(file, list = TRUE)
        for (content in tar_contents) {
          extracted_file <- file.path(dirname(file), content)
          if (file.exists(extracted_file)) {
            extracted_files <- c(extracted_files, extracted_file)
          }
        }
      } else {
        extracted_files <- c(extracted_files, file)
      }
    }
    
    # Use extracted files instead of original downloads
    final_files <- extracted_files
    
    # Check file integrity
    valid_files <- check_file_integrity(final_files)
    if (!all(valid_files)) {
      warning("Some files appear corrupted or empty: ", 
              paste(final_files[!valid_files], collapse = ", "))
      final_files <- final_files[valid_files]
    }
    
    if (verbose) {
      cat("Successfully processed", length(final_files), "files:\n")
      for (f in final_files) {
        cat("  -", basename(f), "\n")
      }
    }
    
    return(final_files)
    
  }, error = function(e) {
    stop("Failed to download GEO data: ", e$message)
  })
}