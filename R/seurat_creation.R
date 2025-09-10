#' Create Seurat Object from Format Information
#'
#' Creates Seurat objects based on detected data format and file organization
#'
#' @param format_info List. Format information from detect_data_format()
#' @param min_cells Numeric. Minimum cells per gene
#' @param min_features Numeric. Minimum features per cell
#' @param project_name Character. Project name for Seurat object
#' @param verbose Logical. Print progress messages
#'
#' @return Seurat object or list of Seurat objects
#' @export
create_seurat_from_format <- function(format_info, 
                                      min_cells = 3,
                                      min_features = 200,
                                      project_name = "scGeoGet",
                                      verbose = TRUE) {
  
  if (is.null(format_info) || is.null(format_info$format)) {
    stop("Invalid format information provided")
  }
  
  format_type <- format_info$format
  
  if (format_type == "10x") {
    return(create_seurat_from_10x(format_info, min_cells, min_features, project_name, verbose))
  }
  
  # Future formats will be added here
  # } else if (format_type == "h5") {
  #   return(create_seurat_from_h5(format_info, min_cells, min_features, project_name, verbose))
  
  stop("Unsupported format type: ", format_type)
}

#' Create Seurat Object from 10X Format
#'
#' @param format_info List with 10X format information
#' @param min_cells Numeric. Minimum cells per gene  
#' @param min_features Numeric. Minimum features per cell
#' @param project_name Character. Project name
#' @param verbose Logical. Print messages
#' @return Seurat object or list of Seurat objects
#' @keywords internal
create_seurat_from_10x <- function(format_info, min_cells, min_features, project_name, verbose) {
  
  samples <- format_info$samples
  
  if (length(samples) == 1) {
    # Single sample case
    if (verbose) cat("Processing single sample...\n")
    
    sample_info <- samples[[1]]
    seurat_obj <- process_10x_sample(sample_info, min_cells, min_features, project_name, verbose)
    return(seurat_obj)
    
  } else {
    # Multiple samples case
    if (verbose) cat("Processing", length(samples), "samples...\n")
    
    seurat_list <- list()
    for (i in seq_along(samples)) {
      sample_name <- names(samples)[i]
      if (verbose) cat("Processing sample:", sample_name, "\n")
      
      sample_project <- paste(project_name, sample_name, sep = "_")
      seurat_list[[sample_name]] <- process_10x_sample(
        samples[[i]], min_cells, min_features, sample_project, verbose
      )
    }
    
    return(seurat_list)
  }
}

#' Process Single 10X Sample
#'
#' @param sample_info List with matrix, barcodes, features file paths
#' @param min_cells Numeric. Minimum cells per gene
#' @param min_features Numeric. Minimum features per cell  
#' @param project_name Character. Project name
#' @param verbose Logical. Print messages
#' @return Seurat object
#' @keywords internal
process_10x_sample <- function(sample_info, min_cells, min_features, project_name, verbose) {
  
  # Validate required files exist
  required_files <- c("matrix", "barcodes", "features")
  missing_files <- setdiff(required_files, names(sample_info))
  if (length(missing_files) > 0) {
    stop("Missing required files for 10X format: ", paste(missing_files, collapse = ", "))
  }
  
  # Check file existence
  file_paths <- unlist(sample_info)
  if (!all(file.exists(file_paths))) {
    missing <- file_paths[!file.exists(file_paths)]
    stop("Files not found: ", paste(missing, collapse = ", "))
  }
  
  tryCatch({
    # Create temporary directory with proper structure for Read10X
    temp_dir <- tempfile()
    dir.create(temp_dir)
    
    # Copy files to temporary directory with standard names
    matrix_dest <- file.path(temp_dir, "matrix.mtx")
    barcodes_dest <- file.path(temp_dir, "barcodes.tsv")
    features_dest <- file.path(temp_dir, "features.tsv")
    
    # Handle compressed files
    copy_and_decompress(sample_info$matrix, matrix_dest)
    copy_and_decompress(sample_info$barcodes, barcodes_dest)
    copy_and_decompress(sample_info$features, features_dest)
    
    if (verbose) cat("Loading count matrix...\n")
    
    # Use Seurat's Read10X function
    count_matrix <- Seurat::Read10X(data.dir = temp_dir)
    
    if (verbose) {
      cat("Matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells\n")
    }
    
    # Create Seurat object
    if (verbose) cat("Creating Seurat object...\n")
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = count_matrix,
      project = project_name,
      min.cells = min_cells,
      min.features = min_features
    )
    
    # Add metadata about source
    seurat_obj@misc$scGeoGet <- list(
      source_files = sample_info,
      creation_date = Sys.time(),
      package_version = utils::packageVersion("scGeoGet")
    )
    
    # Clean up temporary directory
    unlink(temp_dir, recursive = TRUE)
    
    if (verbose) {
      cat("Seurat object created successfully!\n")
      cat("Final dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
    }
    
    return(seurat_obj)
    
  }, error = function(e) {
    # Clean up on error
    if (exists("temp_dir") && dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
    stop("Failed to create Seurat object: ", e$message)
  })
}

#' Copy and Decompress File if Needed
#'
#' @param source Character. Source file path
#' @param dest Character. Destination file path  
#' @keywords internal
copy_and_decompress <- function(source, dest) {
  if (grepl("\\.gz$", source)) {
    # Decompress gzipped file
    con_in <- gzfile(source, "rb")
    con_out <- file(dest, "wb")
    writeBin(readBin(con_in, "raw", file.info(source)$size), con_out)
    close(con_in)
    close(con_out)
  } else {
    # Just copy regular file
    file.copy(source, dest)
  }
}