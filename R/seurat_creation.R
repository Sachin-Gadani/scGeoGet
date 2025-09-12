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
  } else if (format_type == "csv") {
    return(create_seurat_from_csv(format_info, min_cells, min_features, project_name, verbose))
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
    # Keep compression if original files were compressed
    matrix_ext <- if (grepl("\\.gz$", sample_info$matrix)) ".gz" else ""
    barcodes_ext <- if (grepl("\\.gz$", sample_info$barcodes)) ".gz" else ""
    features_ext <- if (grepl("\\.gz$", sample_info$features)) ".gz" else ""
    
    matrix_dest <- file.path(temp_dir, paste0("matrix.mtx", matrix_ext))
    barcodes_dest <- file.path(temp_dir, paste0("barcodes.tsv", barcodes_ext))
    features_dest <- file.path(temp_dir, paste0("features.tsv", features_ext))
    
    # Copy files preserving compression state
    if (verbose) cat("Copying files (preserving compression)...\n")
    
    copy_file_preserve_compression(sample_info$matrix, matrix_dest, verbose = verbose)
    copy_file_preserve_compression(sample_info$barcodes, barcodes_dest, verbose = verbose) 
    copy_file_preserve_compression(sample_info$features, features_dest, verbose = verbose)
    
    # Verify files were created successfully
    if (!file.exists(matrix_dest)) {
      stop("Failed to create matrix file: ", matrix_dest)
    }
    if (!file.exists(barcodes_dest)) {
      stop("Failed to create barcodes file: ", barcodes_dest)
    }
    if (!file.exists(features_dest)) {
      stop("Failed to create features file: ", features_dest)
    }
    
    if (verbose) {
      cat("Files copied successfully:\n")
      cat("  - Matrix:", basename(matrix_dest), "(", file.info(matrix_dest)$size, "bytes )\n")
      cat("  - Barcodes:", basename(barcodes_dest), "(", file.info(barcodes_dest)$size, "bytes )\n") 
      cat("  - Features:", basename(features_dest), "(", file.info(features_dest)$size, "bytes )\n")
    }
    
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

#' Copy File Preserving Compression State
#'
#' @param source Character. Source file path
#' @param dest Character. Destination file path
#' @param verbose Logical. Print progress messages
#' @keywords internal
copy_file_preserve_compression <- function(source, dest, verbose = FALSE) {
  if (!file.exists(source)) {
    stop("Source file does not exist: ", source)
  }
  
  if (verbose) cat("  Copying:", basename(source), "->", basename(dest), "\n")
  
  tryCatch({
    # Simple file copy - preserve compression state
    success <- file.copy(source, dest, overwrite = TRUE)
    if (!success) {
      stop("Failed to copy file from ", source, " to ", dest)
    }
    
    # Verify the destination file was created and has content
    if (!file.exists(dest)) {
      stop("Destination file was not created: ", dest)
    }
    
    dest_size <- file.info(dest)$size
    if (dest_size == 0) {
      stop("Destination file is empty: ", dest)
    }
    
    if (verbose) cat("    Success:", dest_size, "bytes\n")
    
  }, error = function(e) {
    stop("Error copying ", basename(source), ": ", e$message)
  })
}

#' Copy and Decompress File if Needed
#'
#' @param source Character. Source file path
#' @param dest Character. Destination file path
#' @param verbose Logical. Print progress messages
#' @keywords internal
copy_and_decompress <- function(source, dest, verbose = FALSE) {
  if (!file.exists(source)) {
    stop("Source file does not exist: ", source)
  }
  
  if (verbose) cat("  Processing:", basename(source), "->", basename(dest), "\n")
  
  tryCatch({
    if (grepl("\\.gz$", source)) {
      # Decompress gzipped file
      if (verbose) cat("    Decompressing .gz file...\n")
      
      # Use R.utils::gunzip if available, otherwise manual decompression
      if (requireNamespace("R.utils", quietly = TRUE)) {
        R.utils::gunzip(source, destname = dest, overwrite = TRUE, remove = FALSE)
      } else {
        # Manual decompression
        con_in <- gzfile(source, "rb")
        con_out <- file(dest, "wb")
        
        # Read in chunks to handle large files
        chunk_size <- 1024 * 1024  # 1MB chunks
        repeat {
          chunk <- readBin(con_in, "raw", chunk_size)
          if (length(chunk) == 0) break
          writeBin(chunk, con_out)
        }
        
        close(con_in)
        close(con_out)
      }
    } else {
      # Just copy regular file
      if (verbose) cat("    Copying uncompressed file...\n")
      success <- file.copy(source, dest, overwrite = TRUE)
      if (!success) {
        stop("Failed to copy file from ", source, " to ", dest)
      }
    }
    
    # Verify the destination file was created and has content
    if (!file.exists(dest)) {
      stop("Destination file was not created: ", dest)
    }
    
    dest_size <- file.info(dest)$size
    if (dest_size == 0) {
      stop("Destination file is empty: ", dest)
    }
    
    if (verbose) cat("    Success:", dest_size, "bytes\n")
    
  }, error = function(e) {
    stop("Error processing ", basename(source), ": ", e$message)
  })
}

#' Create Seurat Object from CSV Format
#'
#' @param format_info List with CSV format information
#' @param min_cells Numeric. Minimum cells per gene  
#' @param min_features Numeric. Minimum features per cell
#' @param project_name Character. Project name
#' @param verbose Logical. Print messages
#' @return Seurat object or list of Seurat objects
#' @keywords internal
create_seurat_from_csv <- function(format_info, min_cells, min_features, project_name, verbose) {
  
  samples <- format_info$samples
  
  if (length(samples) == 1) {
    # Single sample case
    if (verbose) cat("Processing single CSV sample...\n")
    
    sample_info <- samples[[1]]
    seurat_obj <- process_csv_sample(sample_info, min_cells, min_features, project_name, verbose)
    return(seurat_obj)
    
  } else {
    # Multiple samples case
    if (verbose) cat("Processing", length(samples), "CSV samples...\n")
    
    seurat_list <- list()
    for (i in seq_along(samples)) {
      sample_name <- names(samples)[i]
      if (verbose) cat("Processing sample:", sample_name, "\n")
      
      sample_project <- paste(project_name, sample_name, sep = "_")
      seurat_list[[sample_name]] <- process_csv_sample(
        samples[[i]], min_cells, min_features, sample_project, verbose
      )
    }
    
    return(seurat_list)
  }
}

#' Process Single CSV Sample
#'
#' @param sample_info List with counts file path and optional annotation
#' @param min_cells Numeric. Minimum cells per gene
#' @param min_features Numeric. Minimum features per cell  
#' @param project_name Character. Project name
#' @param verbose Logical. Print messages
#' @return Seurat object
#' @keywords internal
process_csv_sample <- function(sample_info, min_cells, min_features, project_name, verbose) {
  
  # Validate required files exist
  if (!"counts" %in% names(sample_info)) {
    stop("Missing required counts file for CSV format")
  }
  
  counts_file <- sample_info$counts
  if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
  }
  
  tryCatch({
    if (verbose) cat("Loading CSV count matrix...\n")
    
    # Read the count matrix
    if (grepl("\\.gz$", counts_file)) {
      count_matrix <- utils::read.csv(gzfile(counts_file), row.names = 1, check.names = FALSE)
    } else {
      count_matrix <- utils::read.csv(counts_file, row.names = 1, check.names = FALSE)
    }
    
    # Convert to matrix and transpose if needed
    # CSV format typically has genes as rows, cells as columns (which is what Seurat expects)
    count_matrix <- as.matrix(count_matrix)
    
    if (verbose) {
      cat("Matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells\n")
    }
    
    # Convert to sparse matrix for efficiency
    count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
    
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
      package_version = utils::packageVersion("scGeoGet"),
      original_format = "csv"
    )
    
    # TODO: Add annotation data if available
    # if ("annotation" %in% names(sample_info)) {
    #   # Load and add cell metadata
    # }
    
    if (verbose) {
      cat("Seurat object created successfully!\n")
      cat("Final dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
    }
    
    return(seurat_obj)
    
  }, error = function(e) {
    stop("Failed to create Seurat object from CSV: ", e$message)
  })
}