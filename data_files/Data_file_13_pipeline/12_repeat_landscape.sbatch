#!/bin/bash
#SBATCH --job-name=rm_landscape
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/12_rm_landscape_%A.out
# Step 12b: repeat landscape + class table from the two RepeatMasker passes.  usage: sbatch 12_repeat_landscape.sbatch <label>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; L=$1; R=$NV/05_annotation/$L/01_repeats; cd $R
G=$(awk '{s+=$2} END{print s}' $NV/06_release/$L.fa.fai)
python3 $NV/scripts/12_repeat_landscape.py $G repeat_landscape_$L.png repeat_landscape_$L.tsv rm_custom/*.out rm_dfam/*.out --tbl rm_custom/genome.fa.tbl | tee repeat_summary.txt
cp repeat_landscape_$L.png repeat_landscape_$L.tsv repeat_summary.txt $NV/06_release/ 2>/dev/null || true; echo "== DONE $(date)"
