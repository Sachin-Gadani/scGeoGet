# scGeoGet: Single Cell GEO Data Retrieval and Seurat Object Creation

[![R Package](https://img.shields.io/badge/R%20Package-v0.1.0-blue.svg)](https://github.com/yourusername/scGeoGet)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`scGeoGet` is an R package that automatically downloads single cell RNA-seq datasets from GEO (Gene Expression Omnibus) and creates Seurat objects. It handles various data formats and edge cases common in single cell data submissions.

## Current Status: MVP (v0.1.0)

This is the **Minimum Viable Product** version with focused functionality:

✅ **Supported**: 10X format (matrix.mtx, barcodes.tsv, features.tsv/genes.tsv)  
✅ **Supported**: Single sample datasets  
✅ **Supported**: Standard file naming conventions  
✅ **Supported**: Compressed and uncompressed files  

🚧 **Coming Soon**: HDF5 format, multi-sample datasets, custom formats

## Installation

```r
# Install from GitHub (when available)
# devtools::install_github("yourusername/scGeoGet")

# For now, clone and install locally
# git clone https://github.com/yourusername/scGeoGet.git
# install.packages("path/to/scGeoGet", repos = NULL, type = "source")
```

## Dependencies

```r
# Core packages (required)
install.packages(c("GEOquery", "Seurat", "Matrix", "dplyr", "stringr"))

# Optional packages (for future features)
BiocManager::install("rhdf5")
install.packages(c("hdf5r", "R.utils"))
```

## Quick Start

```r
library(scGeoGet)

# Download and create Seurat object from 10X format dataset
seurat_obj <- scGeoGet("GSE123456")

# With custom parameters
seurat_obj <- scGeoGet("GSE123456", 
                      output_dir = "~/data/GSE123456",
                      min_cells = 5,
                      min_features = 100,
                      project_name = "MyProject")
```

## How It Works

1. **Downloads** supplementary files from GEO using `getGEOSuppFiles()`
2. **Detects** data format (currently 10X format only)
3. **Creates** Seurat object using appropriate format-specific functions

## Function Reference

### Main Functions

- `scGeoGet()` - Main interface for downloading and processing GEO datasets
- `detect_data_format()` - Identifies data format from file patterns  
- `download_geo_data()` - Downloads and organizes GEO supplementary files
- `create_seurat_from_format()` - Format-specific Seurat object creation

### Utility Functions

- `validate_geo_accession()` - Validates GEO accession format
- `parse_geo_metadata()` - Extracts metadata from GEO
- `check_file_integrity()` - Validates downloaded files

## Current Limitations (MVP)

- **10X Format Only**: Only supports 10X Genomics format files
- **Single Sample**: Multi-sample datasets may not work optimally
- **Standard Naming**: Expects conventional file naming (matrix.mtx, barcodes.tsv, etc.)
- **Limited Error Recovery**: Basic error handling and fallback strategies

## Development Roadmap

### Phase 1: MVP ✅ (Current)
- 10X format support
- Single sample datasets
- Basic error handling

### Phase 2: Robustness (Next)
- HDF5 format support (`Read10X_h5()`)
- Multi-sample dataset handling
- Improved error messages and recovery

### Phase 3: Edge Cases (Future)
- Custom format support
- Metadata enhancement
- Performance optimization

## Contributing

This package follows a phased development approach. See `scGeoGet_plan.txt` for detailed development guidelines.

When contributing:
1. Check the plan document for current priorities
2. Update the plan document when adding new functionality
3. Test with real GEO datasets
4. Document edge cases and limitations

## Examples

### Basic Usage

```r
# Simple download and processing
pbmc <- scGeoGet("GSE164378")  # 10X PBMC dataset
print(pbmc)
```

### Custom Parameters

```r
# With quality filtering
seurat_obj <- scGeoGet("GSE123456",
                      min_cells = 10,      # Gene must be in ≥10 cells
                      min_features = 500,  # Cell must have ≥500 genes
                      project_name = "MyExperiment")
```

### Error Handling

```r
# The package provides informative error messages
tryCatch({
  seurat_obj <- scGeoGet("GSE999999")  # Non-existent accession
}, error = function(e) {
  cat("Error:", e$message, "\n")
})
```

## Known Issues

- Some datasets may have non-standard file organization
- Multi-sample datasets currently process only the first sample
- Large datasets may require significant memory

## Citation

If you use scGeoGet in your research, please cite:

```
scGeoGet: Single Cell GEO Data Retrieval and Seurat Object Creation
https://github.com/yourusername/scGeoGet
```

## License

MIT License - see LICENSE file for details.

## Support

- 📁 **Issues**: [GitHub Issues](https://github.com/yourusername/scGeoGet/issues)
- 📖 **Documentation**: See function help (`?scGeoGet`)
- 📋 **Development Plan**: See `scGeoGet_plan.txt`