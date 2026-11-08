#!/bin/bash
#SBATCH --job-name=trna_filter
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06d_trna_filter_%A.out
# Step 14c: tRNAscan-SE 2.0 EukHighConfidenceFilter on the raw tRNA set (large vertebrate genomes carry thousands of tRNA-derived SINE copies).
#   usage: sbatch 06d_trna_filter.sbatch <label>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/ncrna; export LC_ALL=C
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; L=$1; cd $NV/05_annotation/$L/04_ncrna
EukHighConfidenceFilter --result trna.out --ss trna.ss --output hc --prefix trna_hc > hc_filter.log 2>&1
cat hc/trna_hc.summary 2>/dev/null | head -n 30
echo "high-confidence tRNA genes: $(grep -vc "^$" hc/trna_hc.out | awk '{print $1-3}')"
awk 'NR>3{n[$5]++} END{for(k in n) print k, n[k]}' hc/trna_hc.out | sort -k2,2nr | head -25 > trna_hc_isotypes.txt; cat trna_hc_isotypes.txt | tr "\n" ";"; echo
echo "== DONE $(date)"
