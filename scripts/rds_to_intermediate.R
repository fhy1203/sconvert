#!/usr/bin/env Rscript
# Convert RDS (Seurat object) to intermediate files
# Uses efficient formats: MTX for sparse, HDF5 for dense

library(Seurat)
library(SeuratObject)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("--input_rds"), type="character", help="Input RDS file path"),
  make_option(c("--output_dir"), type="character", help="Output directory for intermediate files"),
  make_option(c("--spatial"), action="store_true", default=FALSE, help="Enable spatial transcriptomics mode")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opts$input_rds) || is.null(opts$output_dir)) {
  print_help(opt_parser)
  stop("All arguments are required", call.=FALSE)
}

# Create output directory
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# Read RDS file
seurat_obj <- readRDS(opts$input_rds)

cat(paste("Seurat object read successfully!\n"))
cat(paste("  Cells:", ncol(seurat_obj), "\n"))
cat(paste("  Genes:", nrow(seurat_obj), "\n"))

# Check for duplicate gene names
if (anyDuplicated(rownames(seurat_obj))) {
  dup_count <- sum(duplicated(rownames(seurat_obj)))
  cat(paste("\n✗ Error: Gene names are not unique (", dup_count, " duplicates found)\n", sep=""), file=stderr())
  cat("  Please make gene names unique before conversion.\n", file=stderr())
  stop("Duplicate gene names found")
}

# Check for duplicate cell names
if (anyDuplicated(colnames(seurat_obj))) {
  dup_count <- sum(duplicated(colnames(seurat_obj)))
  cat(paste("\n✗ Error: Cell names are not unique (", dup_count, " duplicates found)\n", sep=""), file=stderr())
  cat("  Please make cell names unique before conversion.\n", file=stderr())
  stop("Duplicate cell names found")
}

# Extract expression matrix
available_layers <- Layers(seurat_obj[["RNA"]])
cat(paste("Available layers:", paste(available_layers, collapse=", "), "\n"))

if ("counts" %in% available_layers) {
  matrix_data <- GetAssayData(seurat_obj, layer = "counts")
} else if ("data" %in% available_layers) {
  matrix_data <- GetAssayData(seurat_obj, layer = "data")
} else {
  matrix_data <- GetAssayData(seurat_obj, layer = available_layers[1])
}

cat(paste("Expression matrix:", nrow(matrix_data), "genes x", ncol(matrix_data), "cells\n"))

# Check if sparse
is_sparse <- inherits(matrix_data, c("dgCMatrix", "dgRMatrix", "dgTMatrix", "sparseMatrix"))
cat(paste("Matrix type:", ifelse(is_sparse, "sparse", "dense"), "\n"))

# Save expression matrix based on type
if (is_sparse) {
  # Use MatrixMarket format for sparse matrices
  cat("Saving as MTX format (sparse)...\n")
  Matrix::writeMM(matrix_data, file.path(opts$output_dir, "matrix.mtx"))

  # Save gene and cell names
  write.csv(data.frame(gene_name = rownames(matrix_data)),
            file.path(opts$output_dir, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(cell_name = colnames(matrix_data)),
            file.path(opts$output_dir, "cells.csv"), row.names = FALSE)

  cat("Sparse matrix saved as MTX\n")
} else {
  # Use HDF5 for dense matrices
  cat("Saving as HDF5 format (dense, compressed)...\n")

  if (requireNamespace("hdf5r", quietly = TRUE)) {
    library(hdf5r)
    h5file <- H5File$file.open(file.path(opts$output_dir, "matrix.h5"), mode = "w")
    h5file[["matrix"]] <- matrix_data
    h5file[["genes"]] <- rownames(matrix_data)
    h5file[["cells"]] <- colnames(matrix_data)
    h5file$close_all()
  } else if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)
    h5write(matrix_data, file.path(opts$output_dir, "matrix.h5"), "matrix")
    h5write(rownames(matrix_data), file.path(opts$output_dir, "matrix.h5"), "genes")
    h5write(colnames(matrix_data), file.path(opts$output_dir, "matrix.h5"), "cells")
  } else {
    stop("Neither hdf5r nor rhdf5 package is available. Please install one: install.packages('hdf5r')")
  }

  cat("Dense matrix saved as HDF5\n")
}

# Save matrix type info
writeLines(ifelse(is_sparse, "sparse", "dense"), file.path(opts$output_dir, "matrix_type.txt"))

# Save cell metadata
meta_data <- seurat_obj@meta.data
rownames(meta_data) <- colnames(matrix_data)

if (ncol(meta_data) > 0) {
  write.csv(meta_data, file.path(opts$output_dir, "metadata.csv"))
  cat(paste("Cell metadata saved:", nrow(meta_data), "x", ncol(meta_data), "\n"))
} else {
  meta_data <- data.frame(row.names = colnames(matrix_data))
  write.csv(meta_data, file.path(opts$output_dir, "metadata.csv"))
  cat("Cell metadata has no additional columns, saved cell names only\n")
}

# Save gene information
tryCatch({
  var_data <- seurat_obj[["RNA"]]@meta.features
  if (!is.null(var_data) && nrow(var_data) > 0 && ncol(var_data) > 0) {
    write.csv(var_data, file.path(opts$output_dir, "variables.csv"))
    cat(paste("Gene information saved:", nrow(var_data), "x", ncol(var_data), "\n"))
  } else {
    var_data <- data.frame(row.names = rownames(matrix_data))
    write.csv(var_data, file.path(opts$output_dir, "variables.csv"))
    cat("Gene information has no additional columns, saved gene names only\n")
  }
}, error = function(e) {
  var_data <- data.frame(row.names = rownames(matrix_data))
  write.csv(var_data, file.path(opts$output_dir, "variables.csv"))
  cat("Gene information does not exist, saved gene names only\n")
})

# Handle spatial data if requested
if (opts$spatial) {
  cat("\n[Spatial mode] Processing spatial transcriptomics data...\n")

  # Check if Seurat object has spatial information
  if (length(seurat_obj@images) > 0) {
    # Get the first image/FOV
    library_id <- names(seurat_obj@images)[1]
    fov <- seurat_obj@images[[library_id]]

    # Extract coordinates (handle both FOV objects and lists)
    coords <- NULL
    scalefactors <- NULL

    if (is(fov, "FOV") || is(fov, "SlideSeq") || is(fov, "VisiumV1")) {
      # It's a proper FOV object
      coords <- fov@coordinates
      if (!is.null(fov@scale.factors) && length(fov@scale.factors) > 0) {
        scalefactors <- fov@scale.factors
      }
    } else if (is.list(fov)) {
      # It's a list (our custom format)
      if ("coordinates" %in% names(fov)) {
        coords <- fov$coordinates
      }
      if ("scale.factors" %in% names(fov)) {
        scalefactors <- fov$scale.factors
      }
    }

    if (!is.null(coords)) {
      # Save coordinates (convert column names back to x, y)
      coords_df <- data.frame(x = coords$imagecol, y = coords$imagerow,
                              row.names = rownames(coords))
      write.csv(coords_df, file.path(opts$output_dir, "spatial_coords.csv"))
      cat(paste("  Spatial coordinates saved:", nrow(coords), "spots\n"))

      # Save scale factors if available
      if (!is.null(scalefactors) && length(scalefactors) > 0) {
        jsonlite::write_json(scalefactors, file.path(opts$output_dir, "spatial_scalefactors.json"))
        cat(paste("  Scale factors saved:", paste(names(scalefactors), collapse=", "), "\n"))
      }

      # Save library ID
      writeLines(library_id, file.path(opts$output_dir, "spatial_library_id.txt"))
      cat(paste("  Library ID:", library_id, "\n"))
    } else {
      cat("  Warning: Could not extract coordinates from spatial data\n")
    }
  } else {
    cat("  Warning: --spatial flag set but no spatial information found in Seurat object\n")
  }
}

cat("\nR part completed!\n")
