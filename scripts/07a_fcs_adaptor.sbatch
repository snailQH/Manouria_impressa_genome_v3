#!/bin/bash
#SBATCH --job-name=fcs_adaptor
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=$PROJECT/01.Genome_Assembly/new-version/logs/07a_fcs_adaptor_%A.out
# Step 11a: NCBI FCS-adaptor 0.5.4 (vector/adapter screen, eukaryote database) on one assembly.  usage: sbatch 07a_fcs_adaptor.sbatch <label> <fasta>
set -euo pipefail
NV=$PROJECT/01.Genome_Assembly/new-version; C=$REFDIR/containers
LABEL=$1; FA=$(readlink -f $2); O=$NV/07_contamination/$LABEL/adaptor; mkdir -p $O; cd $O
echo "== $LABEL $(date) node=$(hostname)"
singularity exec -B <shared-storage>,<shared-storage> $C/fcs-adaptor.sif /app/fcs/bin/av_screen_x -o $O --euk $FA > av_screen.log 2>&1
echo "adaptor: $(grep -vc '^#' $O/fcs_adaptor_report.txt || true) flagged sequences"; head -n 20 $O/fcs_adaptor_report.txt
echo "== DONE $(date)"
