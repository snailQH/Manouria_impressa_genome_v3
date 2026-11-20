#!/bin/bash
#SBATCH --job-name=telomere
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=$PROJECT/01.Genome_Assembly/new-version/logs/09_telomere_%A.out
# Step 16a: telomeric repeat arrays at both ends of every chromosome (supports "chromosome-scale" claim).  usage: sbatch 09_telomere.sbatch <label>
set -euo pipefail
source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/assembly
NV=$PROJECT/01.Genome_Assembly/new-version; L=$1; cd $NV/06_release
python3 $NV/scripts/09_telomere.py $L.fa 50000 > $L.telomeres.tsv
awk -F'\t' 'NR>1{n[$7]++} END{for(k in n) print k, n[k]}' $L.telomeres.tsv; echo "== DONE $(date)"
