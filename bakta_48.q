#!/bin/bash
#PBS -N bakta_48
#PBS -P CFP04-CF-013
#PBS -l select=1:ncpus=12:mem=256gb
#PBS -l walltime=48:00:00
#PBS -q auto
#PBS -j oe
#PBS -M eriv141@visitor.nus.edu.sg
#PBS -m abe
#PBS -J 1-105%20

set -euo pipefail

echo "----------------------------------------"
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Array index: ${PBS_ARRAY_INDEX}"
echo "Working directory: $(pwd)"
echo "----------------------------------------"

# Real paths, no symlinks
PROJECT=/scratch/projects/CFP04/CFP04-CF-013/2_Helene
SEQDIR=${PROJECT}/1_data/2_processed/1_metagenomes/3_assembly/1_tpjlc_peat_ssa/
DB=${PROJECT}/3_databases/db-light
WORKDIR=${PROJECT}/working_dir_H/bact_arch_project
OUTDIR=${WORKDIR}/annotated_bact_all

mkdir -p "${WORKDIR}"
mkdir -p "${OUTDIR}"
cd "${WORKDIR}" || { echo "ERROR: Failed to cd into ${WORKDIR}"; exit 1; }

# HOME is fine for conda since 2_miniconda3 is a real directory
source "${HOME}/2_miniconda3/etc/profile.d/conda.sh"
conda activate bakta_env

# Pick sample for this array index
SAMPLE=$(find "${SEQDIR}" -type f -printf '%f\n' \
    | cut -d '_' -f 1,2,3,4 \
    | sort -u \
    | sed -n "${PBS_ARRAY_INDEX}p")

echo "----------------------------------------"
echo "Annotating ${SAMPLE} at $(date)"
echo "Input file: ${SEQDIR}/${SAMPLE}_contigs.fa"
echo "----------------------------------------"

bakta --db "${DB}" \
    --threads 12 \
    --prefix "${SAMPLE}" \
    --output "${OUTDIR}/${SAMPLE}" \
    --meta \
    --skip-sorf \
    --force \
    "${SEQDIR}/${SAMPLE}_contigs.fa"

echo "----------------------------------------"
echo "Job finished at: $(date)"
echo "----------------------------------------"
