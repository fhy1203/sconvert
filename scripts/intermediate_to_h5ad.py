#!/usr/bin/env python3
"""
Convert intermediate files to H5AD (AnnData)
Supports both sparse (MTX) and dense (HDF5) formats
"""
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io as sio
import h5py
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Convert intermediate files to H5AD')
    parser.add_argument('--input_dir', required=True, help='Input directory with intermediate files')
    parser.add_argument('--output_h5ad', required=True, help='Output H5AD file path')
    parser.add_argument('--spatial', action='store_true', help='Enable spatial transcriptomics mode')

    args = parser.parse_args()
    input_dir = Path(args.input_dir)

    try:
        # Read matrix type
        with open(input_dir / "matrix_type.txt", 'r') as f:
            matrix_type = f.read().strip()

        print(f"Matrix type: {matrix_type}")

        # Read expression matrix based on type
        if matrix_type == "sparse":
            # Read sparse matrix from MTX format
            print("Reading sparse matrix from MTX format...")
            matrix_t = sio.mmread(input_dir / "matrix.mtx")

            # Read gene and cell names
            genes = pd.read_csv(input_dir / "genes.csv")['gene_name'].tolist()
            cells = pd.read_csv(input_dir / "cells.csv")['cell_name'].tolist()

            # Convert to CSC format if not already
            if not hasattr(matrix_t, 'tocsr'):
                from scipy.sparse import csc_matrix
                matrix_t = csc_matrix(matrix_t)

            # Transpose to cells x genes (AnnData format)
            X = matrix_t.T.tocsr()

            print(f"Sparse matrix loaded: {len(genes)} genes x {len(cells)} cells")
        else:
            # Read dense matrix from HDF5 format
            print("Reading dense matrix from HDF5 format...")
            with h5py.File(input_dir / "matrix.h5", 'r') as f:
                matrix_t = f['matrix'][:]
                genes = f['genes'][:].astype(str).tolist()
                cells = f['cells'][:].astype(str).tolist()

            # Transpose to cells x genes
            X = matrix_t.T

            print(f"Dense matrix loaded: {len(genes)} genes x {len(cells)} cells")

        # Read cell metadata
        obs = pd.read_csv(input_dir / "metadata.csv", index_col=0)

        # Filter out "Unnamed" columns
        if obs.shape[1] > 0:
            valid_cols = [col for col in obs.columns if not col.startswith('Unnamed')]
            if len(valid_cols) == 0:
                print("Cell metadata has no additional columns, creating empty DataFrame")
                obs = pd.DataFrame(index=cells)
            else:
                obs = obs[valid_cols]
                # Ensure metadata order matches
                if not set(cells).issubset(set(obs.index)):
                    print("Warning: Metadata row names don't match cell names, using cell names from matrix")
                    obs = pd.DataFrame(index=cells)
                else:
                    obs = obs.loc[cells]

        # Read gene information
        var_file = input_dir / "variables.csv"
        if var_file.exists():
            var = pd.read_csv(var_file, index_col=0)

            # Filter out "Unnamed" columns
            if var.shape[1] > 0:
                valid_cols = [col for col in var.columns if not col.startswith('Unnamed')]
                if len(valid_cols) == 0:
                    print("Gene information has no additional columns, creating empty DataFrame")
                    var = pd.DataFrame(index=genes)
                else:
                    var = var[valid_cols]
                    # Ensure gene information order matches
                    if not set(genes).issubset(set(var.index)):
                        print("Warning: Gene information row names don't match gene names, using gene names from matrix")
                        var = pd.DataFrame(index=genes)
                    else:
                        var = var.loc[genes]
        else:
            var = pd.DataFrame(index=genes)

        print(f"Gene information: {var.shape}")

        # Create AnnData object
        adata = ad.AnnData(
            X=X,
            obs=obs,
            var=var
        )

        print(f"\nAnnData object created successfully:")
        print(f"  Cells: {adata.n_obs}")
        print(f"  Genes: {adata.n_vars}")
        print(f"  Matrix type: {'sparse' if hasattr(X, 'toarray') else 'dense'}")
        print(f"  Cell metadata columns: {list(adata.obs.columns)}")
        print(f"  Gene information columns: {list(adata.var.columns)}")

        # Handle spatial data if requested
        if args.spatial:
            print(f"\n[Spatial mode] Restoring spatial transcriptomics data...")

            coords_file = input_dir / "spatial_coords.csv"
            if coords_file.exists():
                # Read spatial coordinates
                coords_df = pd.read_csv(coords_file, index_col=0)

                # Add to obsm
                adata.obsm['spatial'] = coords_df[['x', 'y']].values
                print(f"  Spatial coordinates restored: {adata.obsm['spatial'].shape}")

                # Restore scale factors and library info if available
                scalefactors_file = input_dir / "spatial_scalefactors.json"
                library_id_file = input_dir / "spatial_library_id.txt"

                if scalefactors_file.exists() and library_id_file.exists():
                    import json
                    with open(library_id_file, 'r') as f:
                        library_id = f.read().strip()
                    with open(scalefactors_file, 'r') as f:
                        scalefactors = json.load(f)

                    # Create spatial structure in uns
                    adata.uns['spatial'] = {
                        library_id: {
                            'scalefactors': scalefactors,
                            'images': {}  # We don't save images
                        }
                    }
                    print(f"  Scale factors restored: {list(scalefactors.keys())}")
                    print(f"  Library ID: {library_id}")
            else:
                print("  Warning: --spatial flag set but no spatial coordinates file found")

        # Save as h5ad
        adata.write_h5ad(args.output_h5ad)
        print(f"\nH5AD file saved: {args.output_h5ad}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
