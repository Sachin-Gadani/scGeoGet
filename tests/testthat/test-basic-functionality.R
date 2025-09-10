# Test basic functionality of scGeoGet
library(testthat)
library(scGeoGet)

# Test dataset: GSE111108 (3 Human Lung Adenocarcinoma cell lines, ~4k cells, 14.5MB)
test_accession <- "GSE111108"

test_that("validate_geo_accession works correctly", {
  expect_true(validate_geo_accession("GSE123456"))
  expect_true(validate_geo_accession("GSE111108"))
  expect_false(validate_geo_accession("invalid"))
  expect_false(validate_geo_accession("gse123"))
  expect_false(validate_geo_accession(123))
  expect_false(validate_geo_accession(c("GSE1", "GSE2")))
})

test_that("file integrity checking works", {
  # Create temporary test files
  temp_file1 <- tempfile()
  temp_file2 <- tempfile()
  writeLines("test content", temp_file1)
  # temp_file2 doesn't exist
  
  result <- check_file_integrity(c(temp_file1, temp_file2, "nonexistent"))
  expect_length(result, 3)
  expect_true(result[1])   # temp_file1 exists
  expect_false(result[2])  # temp_file2 doesn't exist
  expect_false(result[3])  # nonexistent doesn't exist
  
  # Clean up
  unlink(temp_file1)
})

# Integration test - requires internet connection
test_that("GSE111108 can be processed (integration test)", {
  skip_on_cran()
  skip_if_offline()
  
  # This is a slower test - only run if explicitly requested
  skip("Integration test - run manually with testthat::test_file()")
  
  temp_dir <- tempdir()
  
  expect_no_error({
    seurat_obj <- scGeoGet(
      accession = test_accession,
      output_dir = file.path(temp_dir, "scGeoGet_test"),
      min_cells = 3,
      min_features = 200,
      verbose = TRUE
    )
  })
  
  expect_s4_class(seurat_obj, "Seurat")
  expect_gt(ncol(seurat_obj), 1000)  # Should have > 1000 cells
  expect_gt(nrow(seurat_obj), 10000) # Should have > 10000 genes
  
  # Check metadata
  expect_true("scGeoGet" %in% names(seurat_obj@misc))
  expect_equal(seurat_obj@misc$scGeoGet$package_version, utils::packageVersion("scGeoGet"))
})