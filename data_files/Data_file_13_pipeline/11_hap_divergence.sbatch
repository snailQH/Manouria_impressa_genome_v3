#!/bin/bash
#SBATCH --job-name=hap_div
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/11_hap_div_%A.out
# Step 16c: divergence between the two haplotype assemblies (and each vs the primary): minimap2 asm5 + paftools call
#   -> SNV / small indel counts and per-Mb densities, aligned fraction. Supports the "two haplotypes" claim with a number.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; O=$NV/04_qc/hap_divergence; mkdir -p $O; cd $O
R=$NV/06_release
for pair in "hap1_vs_hap2 $R/Mimp_v3_hap2.fa $R/Mimp_v3_hap1.fa" "hap1_vs_primary $R/Mimp_v3.fa $R/Mimp_v3_hap1.fa" "hap2_vs_primary $R/Mimp_v3.fa $R/Mimp_v3_hap2.fa"; do
  set -- $pair; N=$1; REF=$2; QRY=$3
  minimap2 -cx asm5 -t 16 --cs $REF $QRY 2>/dev/null | sort -k6,6 -k8,8n > $N.paf
  paftools.js call -l 10000 -L 50000 $N.paf > $N.var 2> $N.call.log
  python3 - $N <<'PY'
import sys; n=sys.argv[1]; snv=ind=cov=0; big=0
for l in open(n+".var"):
    p=l.split("\t")
    if p[0]=="R": cov+=int(p[3])-int(p[2])
    elif p[0]=="V":
        a,b=p[6],p[7]
        if len(a)==1 and len(b)==1 and a!="-" and b!="-": snv+=1
        elif max(len(a),len(b))>=50: big+=1
        else: ind+=1
print(f"{n}\taligned_bp={cov:,}\tSNV={snv:,} ({1000*snv/cov:.2f}/kb)\tsmall_indel={ind:,} ({1000*ind/cov:.3f}/kb)\tSV>=50bp={big:,}")
PY
done | tee hap_divergence_summary.tsv
echo "== DONE $(date)"
