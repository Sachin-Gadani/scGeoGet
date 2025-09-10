#' Validate GEO Accession Format
#'
#' @param accession Character. GEO accession to validate
#' @return Logical. TRUE if valid format, FALSE otherwise
#' @keywords internal
validate_geo_accession <- function(accession) {
  if (!is.character(accession) || length(accession) != 1) {
    return(FALSE)
  }
  
  # Check for GSE format (GSE followed by digits)
  grepl("^GSE\\d+$", accession)
}

#' Parse GEO Metadata
#'
#' Extract basic metadata from GEO accession
#'
#' @param accession Character. GEO accession number
#' @return List containing metadata or NULL if failed
#' @keywords internal
parse_geo_metadata <- function(accession) {
  tryCatch({
    # Try to get basic GEO info (this works even for NGS data)
    gse_info <- GEOquery::getGEO(accession, GSEMatrix = FALSE, getGPL = FALSE)
    
    # Extract basic information
    metadata <- list(
      title = GEOquery::Meta(gse_info)$title,
      summary = GEOquery::Meta(gse_info)$summary,
      organism = GEOquery::Meta(gse_info)$taxon,
      submission_date = GEOquery::Meta(gse_info)$submission_date,
      samples = names(GEOquery::GSMList(gse_info))
    )
    
    return(metadata)
  }, error = function(e) {
    warning("Could not retrieve GEO metadata: ", e$message)
    return(NULL)
  })
}

#' Check File Integrity
#'
#' Basic checks for file existence and readability
#'
#' @param files Character vector of file paths
#' @return Logical vector indicating which files are valid
#' @keywords internal
check_file_integrity <- function(files) {
  if (length(files) == 0) return(logical(0))
  
  sapply(files, function(f) {
    file.exists(f) && file.info(f)$size > 0
  })
}

#' Standardize Gene Names
#'
#' Basic gene name standardization (placeholder for future enhancement)
#'
#' @param genes Character vector of gene names
#' @return Character vector of standardized gene names
#' @keywords internal
standardize_gene_names <- function(genes) {
  # For now, just return as-is
  # Future: convert between different gene ID formats
  genes
}