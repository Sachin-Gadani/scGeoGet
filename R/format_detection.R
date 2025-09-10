#' Detect Data Format from File List
#'
#' Analyzes downloaded files to determine the data format and organize them
#' appropriately. Currently supports 10X format detection (MVP version).
#'
#' @param file_list Character vector of file paths
#'
#' @return List containing format information, or NULL if no supported format detected
#'   \describe{
#'     \item{format}{Character. Format type (currently only "10x")}
#'     \item{samples}{List. Sample information with file mappings}
#'     \item{files}{Named list of file paths organized by type}
#'   }
#' @export
detect_data_format <- function(file_list) {
  
  if (length(file_list) == 0) {
    return(NULL)
  }
  
  # Get just the filenames for pattern matching
  filenames <- basename(file_list)
  
  # Try to detect 10X format
  format_10x <- detect_10x_format(file_list, filenames)
  if (!is.null(format_10x)) {
    return(format_10x)
  }
  
  # Future: Add HDF5 detection, other formats
  # format_h5 <- detect_h5_format(file_list, filenames)
  # if (!is.null(format_h5)) {
  #   return(format_h5)
  # }
  
  return(NULL)
}

#' Detect 10X Format Files
#'
#' @param file_list Character vector of full file paths
#' @param filenames Character vector of basenames
#' @return List with 10X format information or NULL
#' @keywords internal
detect_10x_format <- function(file_list, filenames) {
  
  # 10X format requires: matrix.mtx, barcodes.tsv, features.tsv (or genes.tsv)
  # Files may be compressed (.gz)
  
  # Pattern matching for core 10X files
  matrix_pattern <- "matrix\\.mtx(\\.gz)?$"
  barcodes_pattern <- "barcodes\\.tsv(\\.gz)?$"
  features_pattern <- "features\\.tsv(\\.gz)?$"
  genes_pattern <- "genes\\.tsv(\\.gz)?$"  # older CellRanger versions
  
  # Find files matching patterns
  matrix_files <- file_list[grepl(matrix_pattern, filenames, ignore.case = TRUE)]
  barcodes_files <- file_list[grepl(barcodes_pattern, filenames, ignore.case = TRUE)]
  features_files <- file_list[grepl(features_pattern, filenames, ignore.case = TRUE)]
  genes_files <- file_list[grepl(genes_pattern, filenames, ignore.case = TRUE)]
  
  # Use features.tsv if available, otherwise genes.tsv
  feature_files <- if (length(features_files) > 0) features_files else genes_files
  
  # Check if we have the minimum required files
  if (length(matrix_files) == 0 || length(barcodes_files) == 0 || length(feature_files) == 0) {
    return(NULL)
  }
  
  # For MVP, handle simple case: one set of files (single sample)
  if (length(matrix_files) == 1 && length(barcodes_files) == 1 && length(feature_files) == 1) {
    return(list(
      format = "10x",
      samples = list(
        sample1 = list(
          matrix = matrix_files[1],
          barcodes = barcodes_files[1],
          features = feature_files[1]
        )
      ),
      files = list(
        matrix = matrix_files,
        barcodes = barcodes_files,
        features = feature_files
      )
    ))
  }
  
  # For multiple files, try to group them by sample
  # This is a simplified approach - future versions will be more sophisticated
  if (length(matrix_files) > 1) {
    warning("Multiple sample detection is experimental. Only processing first sample.")
    return(list(
      format = "10x", 
      samples = list(
        sample1 = list(
          matrix = matrix_files[1],
          barcodes = barcodes_files[1], 
          features = feature_files[1]
        )
      ),
      files = list(
        matrix = matrix_files[1],
        barcodes = barcodes_files[1],
        features = feature_files[1]
      )
    ))
  }
  
  return(NULL)
}