#!/bin/bash
# ============================================================
# 01_link_bismark_reports.sh
# Purpose: Create hardlinks to bismark CpG reports for each
#          DMRichR pairwise comparison and copy sample_info.xlsx
#          into each comparison folder
# Input:   results/LPA/10_DMRichR_inputs/{males,females}/sample_ids/
#          results/LPA/10_DMRichR_inputs/{males,females}/sample_info/
# Output:  results/DMRichR/{males,females}/ASD_comparisons/
#          results/DMRichR/{males,females}/ASD_neurotypical_comparisons/
# Usage:   bash scripts/DMRichR/01_link_bismark_reports.sh
# --- Dry run toggle ---
# DRY_RUN=1 to preview without creating links
# DRY_RUN=0 to actually create hardlinks
# DRY_RUN=1

# --- Sex to process ---
# Options: "males" "females" "both"
# SEX="both"

# ============================================================

# --- Project root ---
PROJECT_ROOT="/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George"

# --- Raw bismark CpG reports location ---
BISMARK_RAW_DIR="/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_epigenerator_George/cytosine_reports_combined"

# --- File extension for bismark reports ---
BISMARK_EXT="_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report.txt.gz"

# --- Input: outputs from script 10 ---
INPUT_BASE="${PROJECT_ROOT}/results/LPA/10_DMRichR_inputs"

# --- Output: DMRichR working directories ---
OUTPUT_BASE="${PROJECT_ROOT}/results/DMRichR"

# --- Dry run toggle ---
# DRY_RUN=1 to preview without creating links
# DRY_RUN=0 to actually create hardlinks
DRY_RUN=0

# --- Sex to process ---
# Options: "males" "females" "both"
SEX="both"

# ============================================================
# Helper: classify comparison into subfolder
# ASD vs ASD  → ASD_comparisons
# ASD vs TD   → ASD_neurotypical_comparisons
# ============================================================
get_subfolder() {
  local batch_name=$1
  if [[ "$batch_name" == *"Neurotypical"* ]]; then
    echo "ASD_neurotypical_comparisons"
  else
    echo "ASD_comparisons"
  fi
}

# ============================================================
# Counters
# ============================================================
total_linked=0
total_missing=0
total_comparisons=0

# ============================================================
# Function: process one sex
# ============================================================
process_sex() {
  local sex=$1
  local id_dir="${INPUT_BASE}/${sex}/sample_ids"
  local info_dir="${INPUT_BASE}/${sex}/sample_info"
  local dest_base="${OUTPUT_BASE}/${sex}"

  echo ""
  echo "================================================================"
  echo "Processing: ${sex}"
  echo "Sample ID directory:   ${id_dir}"
  echo "Sample info directory: ${info_dir}"
  echo "Output base:           ${dest_base}"
  echo "================================================================"

  # --- Check directories exist ---
  if [[ ! -d "$id_dir" ]]; then
    echo "ERROR: Sample ID directory not found: ${id_dir}"
    echo "Have you run script 10 (10_prepare_DMRichR_inputs.R) for ${sex}?"
    return 1
  fi

  if [[ ! -d "$info_dir" ]]; then
    echo "ERROR: Sample info directory not found: ${info_dir}"
    echo "Have you run script 10 (10_prepare_DMRichR_inputs.R) for ${sex}?"
    return 1
  fi

  # --- Check there are ID files to process ---
  shopt -s nullglob
  id_files=("${id_dir}"/*_sample_ids.txt)
  shopt -u nullglob

  if [[ ${#id_files[@]} -eq 0 ]]; then
    echo "ERROR: No *_sample_ids.txt files found in: ${id_dir}"
    return 1
  fi

  echo "Found ${#id_files[@]} comparison(s) to process."

  # --- Create subfolders ---
  if [[ $DRY_RUN -eq 0 ]]; then
    mkdir -p "${dest_base}/ASD_comparisons"
    mkdir -p "${dest_base}/ASD_neurotypical_comparisons"
  else
    echo "[DRY RUN] Would create: ${dest_base}/ASD_comparisons"
    echo "[DRY RUN] Would create: ${dest_base}/ASD_neurotypical_comparisons"
  fi

  local sex_linked=0
  local sex_missing=0
  local sex_comparisons=0

  for id_file in "${id_files[@]}"; do
    batch_name=$(basename "$id_file" _sample_ids.txt)

    # --- Determine which subfolder this comparison belongs to ---
    subfolder=$(get_subfolder "$batch_name")
    batch_dir="${dest_base}/${subfolder}/${batch_name}"

    echo ""
    echo "--- Comparison: ${batch_name} ---"
    echo "    Subfolder:  ${subfolder}"

    if [[ $DRY_RUN -eq 1 ]]; then
      echo "[DRY RUN] Would create directory: ${batch_dir}"
    else
      mkdir -p "$batch_dir"
    fi

    # --- Copy sample_info.xlsx into comparison folder ---
    sample_info_src="${info_dir}/${batch_name}_sample_info.xlsx"
    sample_info_dst="${batch_dir}/sample_info.xlsx"

    if [[ -f "$sample_info_src" ]]; then
      if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY RUN] Would copy sample_info.xlsx: ${batch_name}"
      else
        cp "$sample_info_src" "$sample_info_dst"
        echo "📋 Copied sample_info.xlsx: ${batch_name}"
      fi
    else
      echo "⚠️  sample_info.xlsx not found: ${sample_info_src}"
    fi

    # --- Link bismark files ---
    local batch_linked=0
    local batch_missing=0

    while IFS= read -r sample_id; do

      # Skip empty lines
      [[ -z "$sample_id" ]] && continue

      src="${BISMARK_RAW_DIR}/${sample_id}${BISMARK_EXT}"
      dst="${batch_dir}/${sample_id}${BISMARK_EXT}"

      if [[ -f "$src" ]]; then
        if [[ $DRY_RUN -eq 1 ]]; then
          echo "[DRY RUN] Would hardlink: ${sample_id}"
        else
          # Create hardlink — shares inode, no extra disk space
          if ln "$src" "$dst" 2>/dev/null; then
            echo "✅ Linked: ${sample_id}"
            ((batch_linked++))
          else
            # Hardlink fails across filesystems — fall back to symlink
            echo "⚠️  Hardlink failed for ${sample_id} — trying symlink"
            if ln -s "$src" "$dst" 2>/dev/null; then
              echo "✅ Symlinked: ${sample_id}"
              ((batch_linked++))
            else
              echo "❌ Both hardlink and symlink failed: ${sample_id}"
              ((batch_missing++))
            fi
          fi
        fi
      else
        echo "⚠️  File not found: ${src}"
        ((batch_missing++))
      fi

    done < "$id_file"

    echo "  Linked: ${batch_linked} | Missing: ${batch_missing}"
    ((sex_linked      += batch_linked))
    ((sex_missing     += batch_missing))
    ((sex_comparisons++))

  done

  echo ""
  echo "--- ${sex} summary ---"
  echo "  Comparisons processed: ${sex_comparisons}"
  echo "  Files linked:          ${sex_linked}"
  echo "  Files missing:         ${sex_missing}"

  ((total_linked      += sex_linked))
  ((total_missing     += sex_missing))
  ((total_comparisons += sex_comparisons))
}

# ============================================================
# Run
# ============================================================

if [[ $DRY_RUN -eq 1 ]]; then
  echo "================================================================"
  echo "DRY RUN MODE — no files will be linked or copied"
  echo "Set DRY_RUN=0 to create actual hardlinks"
  echo "================================================================"
fi

case "$SEX" in
  "males")
    process_sex "males"
    ;;
  "females")
    process_sex "females"
    ;;
  "both")
    process_sex "males"
    process_sex "females"
    ;;
  *)
    echo "ERROR: SEX must be 'males', 'females', or 'both'"
    echo "Current value: ${SEX}"
    exit 1
    ;;
esac

# ============================================================
# Final summary
# ============================================================
echo ""
echo "================================================================"
echo "All done."
echo "  Total comparisons: ${total_comparisons}"
echo "  Total linked:      ${total_linked}"
echo "  Total missing:     ${total_missing}"
if [[ $DRY_RUN -eq 1 ]]; then
  echo ""
  echo "This was a DRY RUN. Set DRY_RUN=0 to create actual hardlinks."
fi
echo "================================================================"