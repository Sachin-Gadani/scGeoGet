#!/usr/bin/env Rscript

# Test script for GSE111108 dataset
# 3 Human Lung Adenocarcinoma cell lines (H2228, NCI-H1975, HCC827)
# ~4,000 cells, 14.5 MB download, 10X Genomics format

library(scGeoGet)

# Set up test parameters
accession <- "GSE111108"
output_dir <- file.path(getwd(), "test_data", accession)
project_name <- "LungCancer_CellLines"

cat("=== scGeoGet Test: GSE111108 ===\n")
cat("Dataset: 3 Human Lung Adenocarcinoma cell lines\n")
cat("Expected: ~4,000 cells, 10X format\n")
cat("Output directory:", output_dir, "\n\n")

# Test the main function
start_time <- Sys.time()

tryCatch({
  # Run scGeoGet
  seurat_obj <- scGeoGet(
    accession = accession,
    output_dir = output_dir,
    min_cells = 3,           # Gene must be in ≥3 cells
    min_features = 200,      # Cell must have ≥200 genes
    project_name = project_name,
    verbose = TRUE
  )
  
  # Report results
  end_time <- Sys.time()
  processing_time <- end_time - start_time
  
  cat("\n=== SUCCESS! ===\n")
  cat("Processing time:", round(processing_time, 2), attr(processing_time, "units"), "\n")
  cat("Seurat object created successfully!\n\n")
  
  # Print object summary
  cat("=== SEURAT OBJECT SUMMARY ===\n")
  print(seurat_obj)
  
  cat("\n=== DATASET DIMENSIONS ===\n")
  cat("Genes (features):", nrow(seurat_obj), "\n")
  cat("Cells:", ncol(seurat_obj), "\n")
  
  cat("\n=== METADATA INFO ===\n")
  if ("scGeoGet" %in% names(seurat_obj@misc)) {
    cat("Source files:\n")
    print(seurat_obj@misc$scGeoGet$source_files)
    cat("Creation date:", as.character(seurat_obj@misc$scGeoGet$creation_date), "\n")
    cat("Package version:", as.character(seurat_obj@misc$scGeoGet$package_version), "\n")
  }
  
  cat("\n=== BASIC QC METRICS ===\n")
  # Calculate basic metrics
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  cat("Median genes per cell:", median(seurat_obj$nFeature_RNA), "\n")
  cat("Median UMIs per cell:", median(seurat_obj$nCount_RNA), "\n") 
  cat("Median mitochondrial %:", round(median(seurat_obj$percent.mt), 2), "%\n")
  
  # Save results
  results_file <- file.path(output_dir, paste0(project_name, "_seurat.rds"))
  saveRDS(seurat_obj, results_file)
  cat("\nSeurat object saved to:", results_file, "\n")
  
}, error = function(e) {
  cat("\n=== ERROR ===\n")
  cat("Error occurred:", e$message, "\n")
  cat("Check your internet connection and try again.\n")
  
  # Print debug info
  cat("\n=== DEBUG INFO ===\n")
  cat("Working directory:", getwd(), "\n")
  cat("Output directory exists:", dir.exists(dirname(output_dir)), "\n")
  cat("scGeoGet version:", as.character(utils::packageVersion("scGeoGet")), "\n")
  
  # Print traceback
  cat("\nFull error traceback:\n")
  traceback()
})

cat("\n=== Test completed ===\n")