#!/bin/bash
#SBATCH --job-name=Make_mini_data
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=400G
#SBATCH --account=genome-center-grp
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err
#SBATCH --chdir=/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests

set -euo pipefail

###################
# Run Information #
###################
start=$(date +%s)
hostname
date

####################
# Load Environment #
####################
source ~/.profile

##############
# Run Script #
##############

python /quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests/scripts/make_minidata/00_make_comethyl_minidata.py \
  --input_dir /quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_epigenerator_George/cytosine_reports_combined \
  --output_dir /quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests/data/processed/ \
  --sample_list /quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests/data/processed/sample_ids.txt \
  --bed /quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/tests/data/processed/NDBS_regions.bed