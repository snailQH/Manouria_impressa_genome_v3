#!/bin/bash
#SBATCH --job-name=ncrna
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06c_ncrna_%A.out
# Step 14b: non-coding RNA annotation on the release assembly: tRNAscan-SE 2.0.12 (eukaryotic model, -E) and barrnap 0.9 (rRNA, euk).
#   usage: sbatch 06c_ncrna.sbatch <label>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/ncrna
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; L=$1; FA=$NV/06_release/$L.fa; O=$NV/05_annotation/$L/04_ncrna; mkdir -p $O; cd $O
# tRNAscan-SE asks interactively whether to overwrite existing output and, without a terminal, loops forever writing the prompt
# to its log (453 GB on 2026-10-20). Always start from a clean set of output names.
rm -f trna.out trna.ss trna.stats trna.bed trnascan.log
export LC_ALL=C; echo "== $L $(date) node=$(hostname)"; (tRNAscan-SE -h 2>&1 || true) | (grep -m1 -iE "tRNAscan-SE [0-9]" || true); barrnap --version 2>&1 | head -1
tRNAscan-SE -E -H --detail --thread 16 -o trna.out -f trna.ss -m trna.stats -b trna.bed $FA > trnascan.log 2>&1
awk 'NR>3 && $5!="Undet" && $5!="Pseudo"{n++} NR>3{t++} END{print "tRNA genes:",t,"(non-pseudo, determined isotype:",n")"}' trna.out
barrnap --kingdom euk --threads 32 --outseq rrna.fa $FA > rrna.gff 2> barrnap.log
grep -v "^#" rrna.gff | awk -F'\t' '{split($9,a,";"); split(a[1],b,"="); n[b[2]]++} END{for(k in n) print "rRNA",k,n[k]}'
echo "== DONE $(date)"
