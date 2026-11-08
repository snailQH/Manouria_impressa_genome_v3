#!/bin/bash
#SBATCH --job-name=asmqc
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/04_qc_%A.out
# Step 4: assembly QC.  usage: sbatch 04_qc.sbatch <label> <fasta> [<hap2.fa for Merqury trio-free phasing>]
#   assembly-stats; compleasm 0.2.10 (vertebrata_odb12 + sauropsida_odb12, current lineages);
#   BUSCO 5.8 vertebrata_odb10 (same lineage as the v2.4 paper, for direct comparison);
#   Merqury QV / k-mer completeness against the HiFi meryl db.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
LABEL=$1; FA=$(readlink -f $2); HAP2=${3:-}; O=$NV/04_qc/$LABEL; mkdir -p $O; cd $O
conda activate /data1/opt/conda/envs/assembly; assembly-stats -t $FA | tee stats.tsv
python3 $NV/scripts/asm_stats.py $FA | tee stats_full.tsv
conda activate /data1/opt/conda/envs/compleasm
for L in vertebrata_odb12 sauropsida_odb12; do
  compleasm run -a $FA -o compleasm_$L -l $L -t 48 -L /data1/reference/busco_downloads/mb_downloads > compleasm_$L.log 2>&1 || true
  echo "== compleasm $L"; cat compleasm_$L/summary.txt 2>/dev/null
done
conda activate /data1/opt/conda/envs/busco
busco -i $FA -l vertebrata_odb10 -m genome -c 48 -o busco_vertebrata --download_path /data1/reference/busco_downloads/busco_downloads --offline > busco.log 2>&1 || true
grep -E "C:|Complete" busco_vertebrata/short_summary*.txt 2>/dev/null | head -2
conda activate /data1/opt/conda/envs/genomeqc
if [ -n "$HAP2" ]; then merqury.sh $NV/00_data_qc/IT1.k21.meryl $FA $HAP2 merqury > merqury.log 2>&1; else merqury.sh $NV/00_data_qc/IT1.k21.meryl $FA merqury > merqury.log 2>&1; fi
echo "== merqury"; cat merqury.qv merqury.completeness.stats 2>/dev/null
echo "== DONE $(date)"
