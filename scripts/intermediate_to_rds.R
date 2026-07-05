#!/usr/bin/env Rscript
# Convert intermediate files to RDS (Seurat object)
# Supports both sparse (MTX) and dense (HDF5) formats

library(Seurat)
library(SeuratObject)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("--input_dir"), type="character", help="Input directory with intermediate files"),
  make_option(c("--output_rds"), type="character", help="Output RDS file path"),
  make_option(c("--spatial"), action="store_true", default=FALSE, help="Enable spatial transcriptomics mode")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opts$input_dir) || is.null(opts$output_rds)) {
  print_help(opt_parser)
  stop("All arguments are required", call.=FALSE)
}

input_dir <- opts$input_dir

# Read matrix type
matrix_type <- readLines(file.path(input_dir, "matrix_type.txt"))
cat(paste("Matrix type:", matrix_type, "\n"))

# Read expression matrix based on type
if (matrix_type == "sparse") {
  # Read sparse matrix from MTX format
  cat("Reading sparse matrix from MTX format...\n")
  matrix_t <- Matrix::readMM(file.path(input_dir, "matrix.mtx"))

  # Read gene and cell names
  genes <- read.csv(file.path(input_dir, "genes.csv"))$gene_name
  cells <- read.csv(file.path(input_dir, "cells.csv"))$cell_name

  # Set dimnames
  rownames(matrix_t) <- genes
  colnames(matrix_t) <- cells

  cat(paste("Sparse matrix loaded:", nrow(matrix_t), "genes x", ncol(matrix_t), "cells\n"))
} else {
  # Read dense matrix from HDF5 format
  cat("Reading dense matrix from HDF5 format...\n")

  # Check if hdf5r is available, otherwise use rhdf5
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    library(hdf5r)
    h5file <- H5File$file.open(file.path(input_dir, "matrix.h5"), mode = "r")
    matrix_t <- h5file[["matrix"]][]
    genes <- h5file[["genes"]][]
    cells <- h5file[["cells"]][]
    h5file$close_all()

    # Convert from raw to character
    genes <- rawToChar(genes)
    cells <- rawToChar(cells)
  } else if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)
    h5_data <- h5read(file.path(input_dir, "matrix.h5"), "matrix")
    genes <- h5read(file.path(input_dir, "matrix.h5"), "genes")
    cells <- h5read(file.path(input_dir, "matrix.h5"), "cells")

    matrix_t <- h5_data
    genes <- rawToChar(genes)
    cells <- rawToChar(cells)
  } else {
    stop("Neither hdf5r nor rhdf5 package is available. Please install one: install.packages('hdf5r')")
  }

  rownames(matrix_t) <- genes
  colnames(matrix_t) <- cells

  cat(paste("Dense matrix loaded:", nrow(matrix_t), "genes x", ncol(matrix_t), "cells\n"))
}

# Read cell metadata
meta_data <- read.csv(file.path(input_dir, "metadata.csv"), row.names = 1)

# Filter out "Unnamed" columns
if (ncol(meta_data) > 0) {
  valid_cols <- !grepl("^Unnamed", colnames(meta_data))
  meta_data <- meta_data[, valid_cols, drop = FALSE]
}

# Read gene information
var_file <- file.path(input_dir, "variables.csv")
if (file.exists(var_file)) {
  var_data <- read.csv(var_file, row.names = 1)

  # Filter out "Unnamed" columns
  if (ncol(var_data) > 0) {
    valid_cols <- !grepl("^Unnamed", colnames(var_data))
    var_data <- var_data[, valid_cols, drop = FALSE]
  }
  cat(paste("Gene information:", nrow(var_data), "x", ncol(var_data), "\n"))
} else {
  var_data <- data.frame()
}

# Create Seurat object
cat("\nCreating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = matrix_t,
  meta.data = meta_data
)

# Add gene information to Seurat object (only if there are actual columns)
if (ncol(var_data) > 0) {
  tryCatch({
    for (col_name in colnames(var_data)) {
      seurat_obj[["RNA"]][[col_name]] <- var_data[[col_name]]
    }
    cat(paste("Gene information added:", nrow(var_data), "genes,", ncol(var_data), "columns\n"))
  }, error = function(e) {
    cat(paste("Warning: Failed to add gene information:", e$message, "\n"))
  })
} else {
  cat("Gene information has no additional columns, using gene names only\n")
}

# Handle spatial data if requested
if (opts$spatial) {
  cat("\n[Spatial mode] Processing spatial transcriptomics data...\n")

  coords_file <- file.path(input_dir, "spatial_coords.csv")
  if (file.exists(coords_file)) {
    # Read spatial coordinates
    coords_df <- read.csv(coords_file, row.names = 1)

    # Read scale factors if available
    scalefactors_file <- file.path(input_dir, "spatial_scalefactors.json")
    library_id_file <- file.path(input_dir, "spatial_library_id.txt")

    # Prepare coordinates for Seurat
    coords_for_seurat <- coords_df
    colnames(coords_for_seurat) <- c("imagerow", "imagecol")

    if (file.exists(scalefactors_file) && file.exists(library_id_file)) {
      library_id <- trimws(readLines(library_id_file))
      scalefactors <- jsonlite::fromJSON(scalefactors_file)

      # Add spatial information directly (not using FOV class to avoid compatibility issues)
      seurat_obj@images <- list()
      seurat_obj@images[[library_id]] <- list(
        coordinates = coords_for_seurat,
        assay = "RNA",
        key = library_id,
        scale.factors = scalefactors
      )

      cat(paste("  Spatial coordinates added for library:", library_id, "\n"))
      cat(paste("  Scale factors:", paste(names(scalefactors), collapse=", "), "\n"))
    } else {
      # No scale factors, just add coordinates
      seurat_obj@images <- list()
      seurat_obj@images[["sample1"]] <- list(
        coordinates = coords_for_seurat,
        assay = "RNA",
        key = "sample1"
      )

      cat("  Spatial coordinates added (no scale factors found)\n")
    }

    # Verify the data was set
    cat(paste("  Verification - Images slot length:", length(seurat_obj@images), "\n"))
  } else {
    cat("  Warning: --spatial flag set but no spatial coordinates file found\n")
  }
}

cat(paste("\nSeurat object created successfully!\n"))
cat(paste("  Cells:", ncol(seurat_obj), "\n"))
cat(paste("  Genes:", nrow(seurat_obj), "\n"))

# Save as RDS
saveRDS(seurat_obj, opts$output_rds)
cat(paste("RDS file saved:", opts$output_rds, "\n"))
