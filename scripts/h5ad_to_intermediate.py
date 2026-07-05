#!/usr/bin/env python3
"""
Convert H5AD (AnnData) to intermediate files
Uses efficient formats: MTX for sparse, HDF5 for dense
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import h5py
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Convert H5AD to intermediate files')
    parser.add_argument('--input_h5ad', required=True, help='Input H5AD file path')
    parser.add_argument('--output_dir', required=True, help='Output directory for intermediate files')
    parser.add_argument('--layer', help='Layer name to convert (optional, default uses X)', default=None)
    parser.add_argument('--spatial', action='store_true', help='Enable spatial transcriptomics mode')

    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Read h5ad file
        adata = sc.read_h5ad(args.input_h5ad)

        print(f"Read AnnData object:")
        print(f"  Cells: {adata.n_obs}")
        print(f"  Genes: {adata.n_vars}")

        # Check for duplicate gene names
        if not adata.var_names.is_unique:
            dup_count = adata.var_names.duplicated().sum()
            print(f"\n✗ Error: Gene names are not unique ({dup_count} duplicates found)", file=sys.stderr)
            print(f"  Please make gene names unique before conversion.", file=sys.stderr)
            print(f"  You can use: adata.var_names_make_unique()", file=sys.stderr)
            print(f"  Example:", file=sys.stderr)
            print(f"    import scanpy as sc", file=sys.stderr)
            print(f"    adata = sc.read_h5ad('{args.input_h5ad}')", file=sys.stderr)
            print(f"    adata.var_names_make_unique()", file=sys.stderr)
            print(f"    adata.write_h5ad('{args.input_h5ad}')  # Save back", file=sys.stderr)
            sys.exit(1)

        # Check for duplicate cell names
        if not adata.obs_names.is_unique:
            dup_count = adata.obs_names.duplicated().sum()
            print(f"\n✗ Error: Cell names are not unique ({dup_count} duplicates found)", file=sys.stderr)
            print(f"  Please make cell names unique before conversion.", file=sys.stderr)
            sys.exit(1)

        # Extract expression matrix based on user-specified layer
        selected_layer = args.layer

        if selected_layer is None or selected_layer == "X":
            print(f"\nUsing matrix: X (main expression matrix)")
            X = adata.X
        elif selected_layer in adata.layers:
            print(f"\nUsing matrix: layers['{selected_layer}']")
            X = adata.layers[selected_layer]
        else:
            print(f"\nError: Specified layer '{selected_layer}' does not exist", file=sys.stderr)
            print(f"Available layers: {list(adata.layers.keys())}", file=sys.stderr)
            sys.exit(1)

        # Check if sparse
        is_sparse = hasattr(X, 'toarray')
        print(f"\nMatrix type: {'sparse' if is_sparse else 'dense'}")

        # Save expression matrix
        if is_sparse:
            # Use MatrixMarket format for sparse matrices (very efficient)
            matrix_file = output_dir / "matrix.mtx"
            sio.mmwrite(str(matrix_file), X.T)  # Transpose to genes x cells for Seurat
            print(f"  Expression matrix saved as MTX (sparse): {X.shape[1]} genes x {X.shape[0]} cells")

            # Save gene and cell names
            pd.DataFrame({'gene_name': adata.var_names}).to_csv(output_dir / "genes.csv", index=False)
            pd.DataFrame({'cell_name': adata.obs_names}).to_csv(output_dir / "cells.csv", index=False)
        else:
            # Use HDF5 for dense matrices (efficient binary format)
            matrix_file = output_dir / "matrix.h5"
            with h5py.File(matrix_file, 'w') as f:
                # Store as genes x cells for Seurat
                f.create_dataset('matrix', data=X.T, compression='gzip')
                f.create_dataset('genes', data=adata.var_names.astype('S'))
                f.create_dataset('cells', data=adata.obs_names.astype('S'))
            print(f"  Expression matrix saved as HDF5 (dense, compressed): {X.shape[1]} genes x {X.shape[0]} cells")

        # Save cell metadata
        if adata.obs.shape[1] > 0:
            adata.obs.to_csv(output_dir / "metadata.csv")
            print(f"  Cell metadata saved: {adata.obs.shape}")
        else:
            pd.DataFrame(index=adata.obs_names).to_csv(output_dir / "metadata.csv")
            print(f"  Cell metadata saved: cell names only ({len(adata.obs_names)} cells)")

        # Save gene information
        if adata.var.shape[1] > 0:
            adata.var.to_csv(output_dir / "variables.csv")
            print(f"  Gene information saved: {adata.var.shape}")
        else:
            pd.DataFrame(index=adata.var_names).to_csv(output_dir / "variables.csv")
            print(f"  Gene information saved: gene names only ({len(adata.var_names)} genes)")

        # Save matrix type info
        with open(output_dir / "matrix_type.txt", 'w') as f:
            f.write('sparse' if is_sparse else 'dense')

        # Handle spatial data if requested
        if args.spatial:
            print(f"\n[Spatial mode] Processing spatial transcriptomics data...")

            # Check if spatial coordinates exist
            if 'spatial' in adata.obsm:
                coords = adata.obsm['spatial']
                # Save coordinates
                pd.DataFrame(coords, index=adata.obs_names, columns=['x', 'y']).to_csv(output_dir / "spatial_coords.csv")
                print(f"  Spatial coordinates saved: {coords.shape}")

                # Save scale factors if available
                if 'spatial' in adata.uns and len(adata.uns['spatial']) > 0:
                    # Get the first library_id
                    library_id = list(adata.uns['spatial'].keys())[0]
                    spatial_info = adata.uns['spatial'][library_id]

                    if 'scalefactors' in spatial_info:
                        scalefactors = spatial_info['scalefactors']
                        import json
                        with open(output_dir / "spatial_scalefactors.json", 'w') as f:
                            json.dump(scalefactors, f)
                        print(f"  Scale factors saved: {list(scalefactors.keys())}")

                    # Save library_id
                    with open(output_dir / "spatial_library_id.txt", 'w') as f:
                        f.write(library_id)
                    print(f"  Library ID: {library_id}")
            else:
                print("  Warning: --spatial flag set but no spatial coordinates found in adata.obsm['spatial']")

        print("\nPython part completed!")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
