#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
# set -e # Uncomment if you want the script to stop on the first error

# --- Configuration ---
# Get the absolute path to the directory where this script resides
SCRIPT_REAL_PATH=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_REAL_PATH")

# --- DEFINE THE BASE DIRECTORY FOR FIGURES (used for RELATIVE_PATHS) ---
FIGURES_BASE_DIR="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/figures/" # <--- EDIT THIS LINE IF NEEDED

# --- DEFINE YOUR RELATIVE TARGET PATHS HERE ---
# These paths are relative to FIGURES_BASE_DIR
RELATIVE_PATHS=(
  # xenium common genes and samples UMAP
  "embed_panel_restricted_genes_and_samples/10x_5um/NSCLC/lung/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "embed_panel_restricted_genes_and_samples/10x_mm_5um/NSCLC/5k/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"

  # chromium common genes UMAP
  "embed_panel_restricted_genes_scrnaseq/Xenium_NSCLC_5k_lung_chromium_common_genes/matched_combo_standard_lung_specific/umap_RNA_counts_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_Level2.1.png"

  # palette cell types
  "palettes/col_palette_cell_types_combo_Level2.1_legend_vertical.pdf"
)

# --- DEFINE YOUR ABSOLUTE TARGET PATHS HERE ---
# These paths are full, absolute paths and ignore FIGURES_BASE_DIR
ABSOLUTE_PATHS=(
  # xenium 5K vs Xenium Lung mean gene expression, shared genes
  "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/xenium_data_analysis/paper_figures/plots/Figure_2/correlation_5k_lung_mean_gene_expr_matched.pdf"
  
  # xenium 5K vs Chromiun and Xenium Lung vs chromium mean gene expression per cell type, shared genes, sample 0PSV
  "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/xenium_data_analysis/paper_figures/plots/Figure_2/selected_celltype_correlation_mean_gene_expr_5k_lung_chromium_0PSV.pdf"
  # xenium 5K vs Chromiun and Xenium Lung vs chromium mean gene expression per cell type, shared genes, sample 1G73
  "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mbilous/Git/xenium_data_analysis/paper_figures/plots/Figure_2/selected_celltype_correlation_mean_gene_expr_5k_lung_chromium_1G73.pdf"


)

# --- Define Max Filename Length (Common limit, adjust if needed for your filesystem) ---
MAX_FILENAME_LEN=250 # Use 255 for many Linux filesystems, 250 is safer

# --- Function to create links (to avoid code duplication) ---
create_link() {
  local TARGET_PATH="$1"
  local LINK_NAME_METHOD="$2" # 'relative_default', 'relative_fallback', 'absolute_basename'
  local ORIGINAL_INPUT_PATH="$3" # The path string from the array

  local LINK_NAME=""

  # --- Generate Link Name ---
  if [[ "$LINK_NAME_METHOD" == "relative_default" ]]; then
      # Use default method for relative paths (replace '/' with '_')
      LINK_NAME="${ORIGINAL_INPUT_PATH//\//_}"
      LINK_NAME="${LINK_NAME#_}" # Remove potential leading '_'
  elif [[ "$LINK_NAME_METHOD" == "relative_fallback" ]]; then
      # Use fallback for relative paths (basename)
      LINK_NAME=$(basename "$ORIGINAL_INPUT_PATH")
      echo "  Using relative fallback link name (basename): '$LINK_NAME'"
  elif [[ "$LINK_NAME_METHOD" == "absolute_basename" ]]; then
       # Use basename method for absolute paths
      LINK_NAME=$(basename "$TARGET_PATH") # Use basename of the full target path
      echo "  Using absolute link name (basename): '$LINK_NAME'"
  else
      echo "  Internal Error: Invalid link name method '$LINK_NAME_METHOD'. Skipping."
      return 1 # Indicate failure
  fi

  # --- Validate Link Name ---
  if [ -z "$LINK_NAME" ] || [ "$LINK_NAME" == "." ] || [ "$LINK_NAME" == ".." ]; then
      echo "  Error: Could not generate a valid link name for '$ORIGINAL_INPUT_PATH' using method '$LINK_NAME_METHOD'. Skipping."
      return 1
  fi
  if [ ${#LINK_NAME} -gt $MAX_FILENAME_LEN ]; then
      local NAME_SNIPPET="${LINK_NAME:0:50}..."
      echo "  Error: Generated link name ('$NAME_SNIPPET') is too long (${#LINK_NAME} chars > $MAX_FILENAME_LEN). Skipping."
      return 1
  fi

  # Construct the full path for the new symbolic link within the script's directory
  local LINK_PATH="$SCRIPT_DIR/$LINK_NAME"

  # --- Collision Check ---
  if [ -e "$LINK_PATH" ]; then
    if [ -L "$LINK_PATH" ] && [ "$(readlink -f "$LINK_PATH")" == "$(readlink -f "$TARGET_PATH")" ]; then
      # Use readlink -f on both sides to resolve intermediate links and compare canonical paths
      echo "  Skipping: Link '$LINK_NAME' already exists and points correctly."
    else
      echo "  Skipping: A file or different link named '$LINK_NAME' already exists."
      echo "            Existing item path: $LINK_PATH"
      echo "            Current target was: $TARGET_PATH"
    fi
    return 0 # Indicate skipped but not failure
  fi

  # --- Create Symlink ---
  echo "  Creating link: '$LINK_PATH' -> '$TARGET_PATH'"
  ln -s "$TARGET_PATH" "$LINK_PATH"

  # --- Verification ---
  if [ $? -eq 0 ]; then
    echo "  Successfully created symbolic link."
    return 0 # Indicate success
  else
    echo "  Error: Failed to create symbolic link '$LINK_NAME'. Check permissions or filesystem limits."
    echo "         Target was: '$TARGET_PATH'"
    return 1 # Indicate failure
  fi
}


# =============================================
# --- Process RELATIVE Paths ---
# =============================================
echo "=== Processing RELATIVE Paths (relative to $FIGURES_BASE_DIR) ==="
echo "Script directory (where links will be created): $SCRIPT_DIR"
echo "Max link name length: $MAX_FILENAME_LEN"
echo "---"

for RELATIVE_PATH in "${RELATIVE_PATHS[@]}"; do
  echo "Processing relative target: '$RELATIVE_PATH'"

  # Handle edge case: Empty relative path string
  if [ -z "$RELATIVE_PATH" ]; then
      echo "  Skipping: Relative path string is empty."
      echo "---"
      continue
  fi

  # Construct the full absolute path to the target file/directory
  # Remove leading/trailing slashes defensively before joining
  CLEAN_BASE_DIR="${FIGURES_BASE_DIR%/}"
  CLEAN_RELATIVE_PATH="${RELATIVE_PATH#/}"
  FULL_TARGET_PATH="${CLEAN_BASE_DIR}/${CLEAN_RELATIVE_PATH}"

  # Check if the *actual* target file exists using the FULL path
  if [ ! -e "$FULL_TARGET_PATH" ]; then
    echo "  Error: Target path '$FULL_TARGET_PATH' does not exist. Skipping link creation."
    echo "---"
    continue
  fi

  # --- Try Default Link Name Generation (replace '/' with '_') ---
  DEFAULT_LINK_NAME="${RELATIVE_PATH//\//_}"
  DEFAULT_LINK_NAME="${DEFAULT_LINK_NAME#_}" # Remove potential leading '_'

  if [ -n "$DEFAULT_LINK_NAME" ] && [ ${#DEFAULT_LINK_NAME} -le $MAX_FILENAME_LEN ]; then
      # Try creating link with the default name
      create_link "$FULL_TARGET_PATH" "relative_default" "$RELATIVE_PATH"
  else
      # Default name is too long or invalid, try fallback (basename)
      if [ -n "$DEFAULT_LINK_NAME" ]; then # Only print length warning if it was generated
         echo "  Warning: Default name ('${DEFAULT_LINK_NAME:0:50}...') too long (${#DEFAULT_LINK_NAME} chars > $MAX_FILENAME_LEN). Trying fallback."
      fi
      create_link "$FULL_TARGET_PATH" "relative_fallback" "$RELATIVE_PATH"
  fi

  echo "---" # Separator for clarity
done


# =============================================
# --- Process ABSOLUTE Paths ---
# =============================================
echo ""
echo "=== Processing ABSOLUTE Paths ==="
echo "---"

for ABSOLUTE_PATH in "${ABSOLUTE_PATHS[@]}"; do
  echo "Processing absolute target: '$ABSOLUTE_PATH'"

  # Handle edge case: Empty absolute path string
  if [ -z "$ABSOLUTE_PATH" ]; then
      echo "  Skipping: Absolute path string is empty."
      echo "---"
      continue
  fi

  # The target path *is* the absolute path
  TARGET_PATH="$ABSOLUTE_PATH"

  # Check if the target path exists
  if [ ! -e "$TARGET_PATH" ]; then
    echo "  Error: Target path '$TARGET_PATH' does not exist. Skipping link creation."
    echo "---"
    continue
  fi

  # --- Generate Link Name using Basename ONLY for absolute paths ---
  create_link "$TARGET_PATH" "absolute_basename" "$ABSOLUTE_PATH"

  echo "---" # Separator for clarity
done


echo ""
echo "============================================="
echo "Finished processing all paths."
echo "============================================="
exit 0