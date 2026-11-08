#!/bin/bash
#SBATCH --job-name=mito_screen
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/08b_mito_screen_%A.out
# Step 15c: find nuclear-assembly sequences that are (mostly) mitochondrial: align the MitoHiFi mitogenome to primary/hap1/hap2;
#   a sequence with >= 80% of its length covered by mito alignment is listed for removal before deposit (NUMTs inside chromosomes are kept).
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; MT=$NV/07_contamination/Mimp_v3/mito/final_mitogenome.fasta; cd $NV/06_release
[ -s $MT ] || { echo "no mitogenome yet"; exit 1; }
# release copy of the mitogenome: FASTA renamed, GenBank annotation, MitoHiFi stats + circular map (Data Records object #6)
M=$(dirname $MT); seqkit replace -p '.*' -r 'Mimp_v3_mito' $MT > Mimp_v3_mito.fa; cp $M/final_mitogenome.gb Mimp_v3_mito.gb 2>/dev/null || true
cp $M/final_mitogenome.annotation.png Mimp_v3_mito.annotation.png 2>/dev/null || true; cp $M/contigs_stats.tsv Mimp_v3_mito.contigs_stats.tsv 2>/dev/null || true
seqkit stats Mimp_v3_mito.fa | tee Mimp_v3_mito.stats.txt; grep -c "gene" Mimp_v3_mito.gb 2>/dev/null | sed "s/^/gene features in GenBank: /"
for A in Mimp_v3 Mimp_v3_hap1 Mimp_v3_hap2; do
  minimap2 -x asm10 -t 16 $A.fa $MT 2>/dev/null | awk '$12>=30' > $A.mito.paf
  python3 - $A <<'PY'
import sys,collections; a=sys.argv[1]; cov=collections.defaultdict(set); L={}
for l in open(a+".mito.paf"):
    p=l.split("\t"); L[p[5]]=int(p[6])
    for i in range(int(p[7])//1000,int(p[8])//1000+1): cov[p[5]].add(i)
with open(a+".mito_hits.tsv","w") as f:
    f.write("sequence\tlength\tmito_covered_kb\tpct\tverdict\n")
    for s in sorted(cov,key=lambda x:-len(cov[x])):
        pct=100*len(cov[s])*1000/L[s]; f.write(f"{s}\t{L[s]}\t{len(cov[s])}\t{pct:.1f}\t{'REMOVE (mitochondrial)' if pct>=80 else 'keep (NUMT)'}\n")
print(a, ": sequences with mito hits", len(cov), "; to remove:", sum(1 for s in cov if 100*len(cov[s])*1000/L[s]>=80))
PY
done; echo "== DONE $(date)"
