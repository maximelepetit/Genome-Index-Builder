#!/bin/bash
set -e

# Script: run_full_mapping_pipeline.sh
# Description: Execute the cell type mapping pipeline from an unlabeled h5ad file.
# Author: 
# Requirements: Python environment with necessary dependencies installed.

# -------------------------------------------------------------------
# Help message
# -------------------------------------------------------------------
show_help() {
  echo "Usage: $0 --input_file INPUT_FILE --output_dir OUTPUT_DIR [options]"
  echo ""
  echo "Required parameters:"
  echo "  --input_file                      Path to the input h5ad file (unlabeled)."
  echo "  --output_dir                      Directory to save the output files."
  echo ""
  echo "Optional parameters:"
  echo "  --n_procs                         Number of processes to use (default=4)."

  exit 1
}

# -------------------------------------------------------------------
# Default values
# -------------------------------------------------------------------
N_PROCS=15

# -------------------------------------------------------------------
# Parse arguments
# -------------------------------------------------------------------
while [[ -n "$1" ]]; do
    case "$1" in
        --input_file) INPUT_FILE="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --n_procs) N_PROCS="$2"; shift 2 ;;
        *) echo "$(date) Unknown option: $1"; show_help ;;
    esac
done

# -------------------------------------------------------------------
# Validate required parameters
# -------------------------------------------------------------------
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" ]]; then
  echo "Error: Missing required parameters."
  echo "You must provide --input_file INPUT_FILE and --output_dir OUTPUT_DIR."
  show_help
fi

# -------------------------------------------------------------------
# Environment setup
# -------------------------------------------------------------------
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"

# -------------------------------------------------------------------
# Preprocessing steps (1–3) if requested
# -------------------------------------------------------------------
echo "Running preprocessing steps..."

  # Step 1: Precompute statistics
  python -m cell_type_mapper.cli.precompute_stats_scrattch \
    --hierarchy '["C7_named", "C25_named", "C66_named", "C185_named", "C286_named"]' \
    --h5ad_path "$INPUT_FILE" \
    --n_processors "$N_PROCS" \
    --normalization raw \
    --clobber True \
    --output_path "$OUTPUT_DIR/precomputed_stats.h5" \
    || { echo "Error: Precompute statistics failed."; exit 1; }

  # Step 2: Select reference markers
  python -m cell_type_mapper.cli.reference_markers \
    --precomputed_path_list "['$OUTPUT_DIR/precomputed_stats.h5']" \
    --n_valid 20 \
    --n_processors "$N_PROCS" \
    --output_dir "$OUTPUT_DIR" \
    --clobber True \
    || { echo "Error: Select reference markers failed."; exit 1; }

  # Step 3: Select query markers
  python -m cell_type_mapper.cli.query_markers \
    --output_path "$OUTPUT_DIR/query_markers.json" \
    --reference_marker_path_list "['$OUTPUT_DIR/reference_markers.h5']" \
    --n_per_utility 10 \
    --n_processors "$N_PROCS" \
    || { echo "Error: Select query markers failed."; exit 1; }

  # Use outputs from preprocessing
  PRECOMPUTED_STATS_PATH="$OUTPUT_DIR/precomputed_stats.h5"
  QUERY_MARKERS_LOOKUP="$OUTPUT_DIR/query_markers.json"

