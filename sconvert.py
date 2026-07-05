#!/usr/bin/env python3
"""
sconvert - Single-cell data format conversion tool
Convert between RDS (Seurat) and H5AD (AnnData) formats
"""

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path

# Get script directory
SCRIPT_DIR = Path(__file__).parent.absolute()
SCRIPTS_DIR = SCRIPT_DIR / "scripts"


def get_rscript_path(custom_path=None):
    """Get Rscript path

    Args:
        custom_path: User-specified Rscript path (optional)

    Returns:
        Rscript path
    """
    # If user specified a path, use it directly
    if custom_path:
        if os.path.exists(custom_path):
            return custom_path
        else:
            raise FileNotFoundError(f"Specified Rscript path does not exist: {custom_path}")

    # Try to use Rscript from current Python environment
    python_path = Path(sys.executable)
    rscript_path = python_path.parent / "Rscript"
    if rscript_path.exists():
        return str(rscript_path)

    # If not found, try to find in PATH
    return "Rscript"


def check_file_exists(filepath):
    """Check if file exists"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File does not exist: {filepath}")


def list_h5ad_layers(h5ad_file):
    """
    List all layers in H5AD file

    Args:
        h5ad_file: H5AD file path

    Returns:
        layers list and X matrix info
    """
    import scanpy as sc

    adata = sc.read_h5ad(h5ad_file)

    layers_info = {
        'X': {
            'shape': adata.X.shape if adata.X is not None else None
        },
        'layers': {}
    }

    # Check layers
    if hasattr(adata, 'layers') and adata.layers:
        for layer_name in adata.layers.keys():
            layer_data = adata.layers[layer_name]
            layers_info['layers'][layer_name] = {
                'shape': layer_data.shape
            }

    return layers_info


def h5ad_to_rds(h5ad_file, rds_file, layer=None, temp_dir=None, rscript_path=None, spatial=False):
    """
    Convert H5AD (AnnData) to RDS (Seurat)

    Args:
        h5ad_file: Input H5AD file path
        rds_file: Output RDS file path
        layer: Layer name to convert, if None use adata.X
        temp_dir: Temporary file directory
        rscript_path: Rscript path (optional)
        spatial: Enable spatial transcriptomics mode (default: False)

    Workflow:
    1. Python reads h5ad -> exports intermediate files (MTX/HDF5)
    2. R reads intermediate files -> creates Seurat object -> saves as rds
    """
    check_file_exists(h5ad_file)

    # Create temporary directory
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="sconvert_")
        cleanup = True
    else:
        os.makedirs(temp_dir, exist_ok=True)
        cleanup = False

    try:
        # Check and display layers in H5AD file
        print(f"[INFO] Analyzing expression matrix in H5AD file...")
        layers_info = list_h5ad_layers(h5ad_file)

        print(f"\nAvailable expression matrices:")
        print(f"  X (main matrix): {layers_info['X']['shape']}")

        if layers_info['layers']:
            for layer_name, info in layers_info['layers'].items():
                print(f"  layers['{layer_name}']: {info['shape']}")
        else:
            print(f"  (no other layers)")

        # Determine which matrix to use
        if layer is None:
            selected_layer = "X"
            print(f"\n[TIP] It is recommended to convert raw count matrix")
            print(f"[TIP] Normalization and log1p processing may differ between R and Python")
            print(f"[INFO] --layer not specified, using default main matrix X")
        else:
            selected_layer = layer
            if layer in layers_info['layers']:
                print(f"\n[INFO] Using specified layer: '{layer}'")
            else:
                print(f"\n[WARNING] Specified layer '{layer}' does not exist, using main matrix X")
                selected_layer = "X"

        print(f"\n[INFO] Converting matrix: {'X (main expression matrix)' if selected_layer == 'X' else f'layers[\'{selected_layer}\']'}")
        print(f"[INFO] Note: Only converting raw expression matrix, no additional data processing")

        print(f"\n[Step 1/2] Python: Reading {h5ad_file} and exporting intermediate files...")

        # Execute Python script to convert H5AD to intermediate files
        py_script = str(SCRIPTS_DIR / "h5ad_to_intermediate.py")
        py_cmd = [
            sys.executable, py_script,
            "--input_h5ad", h5ad_file,
            "--output_dir", temp_dir
        ]

        if selected_layer != "X":
            py_cmd.extend(["--layer", selected_layer])

        if spatial:
            py_cmd.append("--spatial")

        result = subprocess.run(py_cmd, capture_output=True, text=True)

        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        if result.returncode != 0:
            raise RuntimeError(f"Python conversion failed: {result.stderr}")

        print(f"\n[Step 2/2] R: Reading intermediate files and creating Seurat object...")

        # Execute R script to convert intermediate files to RDS
        rscript_exec = get_rscript_path(rscript_path)
        r_script = str(SCRIPTS_DIR / "intermediate_to_rds.R")
        r_cmd = [
            rscript_exec, r_script,
            "--input_dir", temp_dir,
            "--output_rds", rds_file
        ]

        if spatial:
            r_cmd.append("--spatial")

        r_result = subprocess.run(r_cmd, capture_output=True, text=True)

        print(r_result.stdout)
        if r_result.stderr:
            print(r_result.stderr, file=sys.stderr)

        if r_result.returncode != 0:
            raise RuntimeError(f"R conversion failed: {r_result.stderr}")

        print(f"\n✓ Conversion completed: {h5ad_file} -> {rds_file}")
        return True

    finally:
        if cleanup and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def rds_to_h5ad(rds_file, h5ad_file, temp_dir=None, rscript_path=None, spatial=False):
    """
    Convert RDS (Seurat) to H5AD (AnnData)

    Args:
        rds_file: Input RDS file path
        h5ad_file: Output H5AD file path
        temp_dir: Temporary file directory
        rscript_path: Rscript path (optional)
        spatial: Enable spatial transcriptomics mode (default: False)

    Workflow:
    1. R reads rds -> exports intermediate files (MTX/HDF5)
    2. Python reads intermediate files -> creates AnnData object -> saves as h5ad
    """
    check_file_exists(rds_file)

    # Create temporary directory
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="sconvert_")
        cleanup = True
    else:
        os.makedirs(temp_dir, exist_ok=True)
        cleanup = False

    try:
        print(f"[Step 1/2] R: Reading {rds_file} and exporting intermediate files...")

        # Execute R script to convert RDS to intermediate files
        rscript_exec = get_rscript_path(rscript_path)
        r_script = str(SCRIPTS_DIR / "rds_to_intermediate.R")
        r_cmd = [
            rscript_exec, r_script,
            "--input_rds", rds_file,
            "--output_dir", temp_dir
        ]

        if spatial:
            r_cmd.append("--spatial")

        r_result = subprocess.run(r_cmd, capture_output=True, text=True)

        print(r_result.stdout)
        if r_result.stderr:
            print(r_result.stderr, file=sys.stderr)

        if r_result.returncode != 0:
            raise RuntimeError(f"R conversion failed: {r_result.stderr}")

        print(f"\n[Step 2/2] Python: Reading intermediate files and creating AnnData object...")

        # Execute Python script to convert intermediate files to H5AD
        py_script = str(SCRIPTS_DIR / "intermediate_to_h5ad.py")
        py_cmd = [
            sys.executable, py_script,
            "--input_dir", temp_dir,
            "--output_h5ad", h5ad_file
        ]

        if spatial:
            py_cmd.append("--spatial")

        result = subprocess.run(py_cmd, capture_output=True, text=True)

        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        if result.returncode != 0:
            raise RuntimeError(f"Python conversion failed: {result.stderr}")

        print(f"\n✓ Conversion completed: {rds_file} -> {h5ad_file}")
        return True

    finally:
        if cleanup and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="sconvert - Single-cell data format conversion tool (RDS <-> H5AD)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # H5AD to RDS (using default X matrix)
  python sconvert.py h5ad2rds --input_h5ad input.h5ad --output_rds output.rds --rscript_path /path/to/Rscript

  # H5AD to RDS (specify layer)
  python sconvert.py h5ad2rds --input_h5ad input.h5ad --output_rds output.rds --rscript_path /path/to/Rscript --layer counts

  # RDS to H5AD
  python sconvert.py rds2h5ad --input_rds input.rds --output_h5ad output.h5ad --rscript_path /path/to/Rscript

  # Specify temp directory
  python sconvert.py h5ad2rds --input_h5ad input.h5ad --output_rds output.rds --rscript_path /path/to/Rscript --temp_dir /tmp/work
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="Conversion command")

    # H5AD to RDS subcommand
    parser_h5ad2rds = subparsers.add_parser("h5ad2rds", help="Convert H5AD (AnnData) to RDS (Seurat)")
    parser_h5ad2rds.add_argument("--input_h5ad", required=True, help="Input H5AD file path")
    parser_h5ad2rds.add_argument("--output_rds", required=True, help="Output RDS file path")
    parser_h5ad2rds.add_argument("--rscript_path", required=True, help="Path to Rscript executable")
    parser_h5ad2rds.add_argument("--layer", help="Specify layer name to convert (optional, default uses X main matrix)", default=None)
    parser_h5ad2rds.add_argument("--temp_dir", help="Temporary file directory (optional)", default=None)
    parser_h5ad2rds.add_argument("--spatial", action="store_true", help="Enable spatial transcriptomics mode (preserve coordinates)")

    # RDS to H5AD subcommand
    parser_rds2h5ad = subparsers.add_parser("rds2h5ad", help="Convert RDS (Seurat) to H5AD (AnnData)")
    parser_rds2h5ad.add_argument("--input_rds", required=True, help="Input RDS file path")
    parser_rds2h5ad.add_argument("--output_h5ad", required=True, help="Output H5AD file path")
    parser_rds2h5ad.add_argument("--rscript_path", required=True, help="Path to Rscript executable")
    parser_rds2h5ad.add_argument("--temp_dir", help="Temporary file directory (optional)", default=None)
    parser_rds2h5ad.add_argument("--spatial", action="store_true", help="Enable spatial transcriptomics mode (preserve coordinates)")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Get rscript_path from subcommand
    rscript_path = args.rscript_path

    try:
        if args.command == "h5ad2rds":
            h5ad_to_rds(
                args.input_h5ad,
                args.output_rds,
                layer=args.layer,
                temp_dir=args.temp_dir,
                rscript_path=rscript_path,
                spatial=args.spatial
            )
        elif args.command == "rds2h5ad":
            rds_to_h5ad(
                args.input_rds,
                args.output_h5ad,
                temp_dir=args.temp_dir,
                rscript_path=rscript_path,
                spatial=args.spatial
            )

    except Exception as e:
        print(f"\n✗ Conversion failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
