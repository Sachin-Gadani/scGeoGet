# Test format detection functionality
library(testthat)
library(scGeoGet)

test_that("10X format detection works correctly", {
  # Test with typical 10X filenames
  files_10x_v3 <- c(
    "/path/to/matrix.mtx.gz",
    "/path/to/barcodes.tsv.gz", 
    "/path/to/features.tsv.gz"
  )
  
  files_10x_v2 <- c(
    "/path/to/matrix.mtx",
    "/path/to/barcodes.tsv",
    "/path/to/genes.tsv"
  )
  
  # Mock the files exist for testing
  old_file_exists <- file.exists
  assign("file.exists", function(x) rep(TRUE, length(x)), envir = .GlobalEnv)
  
  # Test v3 format detection
  result_v3 <- detect_data_format(files_10x_v3)
  expect_equal(result_v3$format, "10x")
  expect_length(result_v3$samples, 1)
  expect_equal(names(result_v3$samples), "sample1")
  
  # Test v2 format detection  
  result_v2 <- detect_data_format(files_10x_v2)
  expect_equal(result_v2$format, "10x")
  expect_length(result_v2$samples, 1)
  
  # Test with non-10X files
  files_other <- c("/path/to/data.csv", "/path/to/metadata.txt")
  result_other <- detect_data_format(files_other)
  expect_null(result_other)
  
  # Test with empty file list
  result_empty <- detect_data_format(c())
  expect_null(result_empty)
  
  # Restore original file.exists
  assign("file.exists", old_file_exists, envir = .GlobalEnv)
})

test_that("10X format detection handles missing files", {
  # Test incomplete 10X format (missing barcodes)
  files_incomplete <- c(
    "/path/to/matrix.mtx.gz",
    "/path/to/features.tsv.gz"
    # missing barcodes.tsv.gz
  )
  
  result <- detect_data_format(files_incomplete)
  expect_null(result)
})

test_that("10X format detection handles mixed compression", {
  # Test mixed compression states
  files_mixed <- c(
    "/path/to/matrix.mtx.gz",      # compressed
    "/path/to/barcodes.tsv",       # uncompressed  
    "/path/to/features.tsv.gz"     # compressed
  )
  
  # Mock file.exists
  old_file_exists <- file.exists
  assign("file.exists", function(x) rep(TRUE, length(x)), envir = .GlobalEnv)
  
  result <- detect_data_format(files_mixed)
  expect_equal(result$format, "10x")
  expect_length(result$samples, 1)
  
  # Restore original file.exists
  assign("file.exists", old_file_exists, envir = .GlobalEnv)
})