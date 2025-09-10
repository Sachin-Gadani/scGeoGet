#' Download Single Cell Dataset from GEO and Create Seurat Object
#'
#' Main interface function for downloading single cell RNA-seq datasets from GEO
#' and creating Seurat objects. Currently supports 10X format (MVP version).
#'
#' @param accession Character. GEO accession number (e.g., "GSE123456")
#' @param output_dir Character. Directory to store downloaded files. If NULL, 
#'   creates a temporary directory.
#' @param min_cells Numeric. Minimum number of cells expressing a gene (default: 3)
#' @param min_features Numeric. Minimum number of features per cell (default: 200)
#' @param project_name Character. Project name for Seurat object (default: accession)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A Seurat object or list of Seurat objects (for multi-sample datasets)
#' 
#' @details
#' This function follows a three-step process:
#' \enumerate{
#'   \item Downloads supplementary files from GEO using \code{getGEOSuppFiles()}
#'   \item Detects data format and organizes files
#'   \item Creates Seurat object(s) using appropriate format-specific functions
#' }
#' 
#' Current limitations (MVP version):
#' \itemize{
#'   \item Supports 10X format (matrix.mtx, barcodes.tsv, features.tsv/genes.tsv)
#'   \item Supports CSV format (gene_count.csv files)
#'   \item Single sample datasets work best
#'   \item Standard file naming conventions expected
#' }
#'
#' @examples
#' \dontrun{
#' # Download and process a 10X format dataset
#' seurat_obj <- scGeoGet("GSE123456")
#' 
#' # With custom parameters
#' seurat_obj <- scGeoGet("GSE123456", 
#'                       output_dir = "~/data/GSE123456",
#'                       min_cells = 5,
#'                       min_features = 100)
#' }
#'
#' @export
scGeoGet <- function(accession, 
                     output_dir = NULL,
                     min_cells = 3,
                     min_features = 200,
                     project_name = NULL,
                     verbose = TRUE) {
  
  # Input validation
  if (!validate_geo_accession(accession)) {
    stop("Invalid GEO accession format. Expected format: GSE######")
  }
  
  if (verbose) cat("Starting scGeoGet for accession:", accession, "\n")
  
  # Set up output directory
  if (is.null(output_dir)) {
    output_dir <- tempdir()
    if (verbose) cat("Using temporary directory:", output_dir, "\n")
  }
  
  # Set default project name
  if (is.null(project_name)) {
    project_name <- accession
  }
  
  tryCatch({
    # Step 1: Download GEO data
    if (verbose) cat("Step 1: Downloading GEO supplementary files...\n")
    downloaded_files <- download_geo_data(accession, output_dir, verbose = verbose)
    
    # Step 2: Detect data format
    if (verbose) cat("Step 2: Detecting data format...\n")
    format_info <- detect_data_format(downloaded_files)
    
    if (is.null(format_info)) {
      stop("Could not detect supported data format. Currently only 10X format is supported.")
    }
    
    # Step 3: Create Seurat object(s)
    if (verbose) cat("Step 3: Creating Seurat object(s)...\n")
    seurat_result <- create_seurat_from_format(
      format_info = format_info,
      min_cells = min_cells,
      min_features = min_features, 
      project_name = project_name,
      verbose = verbose
    )
    
    if (verbose) cat("Successfully created Seurat object(s)!\n")
    return(seurat_result)
    
  }, error = function(e) {
    stop("Error in scGeoGet: ", e$message, call. = FALSE)
  })
}