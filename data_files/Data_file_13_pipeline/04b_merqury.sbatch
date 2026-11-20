#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=$PROJECT/01.Genome_Assembly/new-version/logs/04b_merqury_%A.out
# Merqury only (k=21 HiFi db from step 1).  usage: sbatch 04b_merqury.sbatch <label> <asm1.fa> [<asm2.fa>]
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/genomeqc
NV=$PROJECT/01.Genome_Assembly/new-version
LABEL=$1; A1=$(readlink -f $2); A2=${3:+$(readlink -f $3)}; O=$NV/04_qc/$LABEL/merqury_run; rm -rf $O; mkdir -p $O; cd $O
export MERQURY=$CONDA_ROOT/envs/genomeqc/share/merqury
ln -sf $NV/00_data_qc/IT1.k21.meryl reads.k21.meryl; ln -sf $A1 asm1.fa; [ -n "$A2" ] && ln -sf $A2 asm2.fa
merqury.sh reads.k21.meryl asm1.fa ${A2:+asm2.fa} merqury > merqury.log 2>&1
echo "== QV"; cat merqury.qv; echo "== completeness"; cat merqury.completeness.stats; echo "== DONE $(date)"
