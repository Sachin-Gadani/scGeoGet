#' scGeoGet: Single Cell GEO Data Retrieval and Seurat Object Creation
#'
#' Downloads single cell RNA-seq datasets from GEO (Gene Expression Omnibus) 
#' and creates Seurat objects. Handles various data formats and edge cases 
#' common in single cell data submissions.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{scGeoGet}}}{Main interface function for downloading and processing GEO datasets}
#'   \item{\code{\link{detect_data_format}}}{Identifies data format from file patterns}
#'   \item{\code{\link{download_geo_data}}}{Downloads and organizes GEO supplementary files}
#'   \item{\code{\link{create_seurat_from_format}}}{Format-specific Seurat object creation}
#' }
#'
#' @section Development Approach:
#' This package follows a phased development approach:
#' \enumerate{
#'   \item MVP: 10X format support only
#'   \item Robustness: Add HDF5 and multi-sample support  
#'   \item Edge Cases: Community-driven expansion
#' }
#'
#' @docType package
#' @name scGeoGet-package
#' @aliases scGeoGet-package
NULL

# Global variables to avoid R CMD check notes
utils::globalVariables(c(".", "file_type", "sample_id"))