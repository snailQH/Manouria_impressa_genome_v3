#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/04b_merqury_%A.out
# Merqury only (k=21 HiFi db from step 1).  usage: sbatch 04b_merqury.sbatch <label> <asm1.fa> [<asm2.fa>]
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/genomeqc
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
LABEL=$1; A1=$(readlink -f $2); A2=${3:+$(readlink -f $3)}; O=$NV/04_qc/$LABEL/merqury_run; rm -rf $O; mkdir -p $O; cd $O
export MERQURY=/data1/opt/conda/envs/genomeqc/share/merqury
ln -sf $NV/00_data_qc/IT1.k21.meryl reads.k21.meryl; ln -sf $A1 asm1.fa; [ -n "$A2" ] && ln -sf $A2 asm2.fa
merqury.sh reads.k21.meryl asm1.fa ${A2:+asm2.fa} merqury > merqury.log 2>&1
echo "== QV"; cat merqury.qv; echo "== completeness"; cat merqury.completeness.stats; echo "== DONE $(date)"
