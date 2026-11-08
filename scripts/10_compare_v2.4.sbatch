#!/bin/bash
#SBATCH --job-name=cmp_v24
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/10_compare_v24_%A.out
# Step 16b: v2.4 (GWHHIVU00000000.1) vs v3 correspondence: minimap2 asm5 of v2.4 sequences onto Mimp_v3 -> per-sequence best match
#   (liftover guide for users of v2.4; Usage Notes + Supplementary Table).  usage: sbatch 10_compare_v2.4.sbatch
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; O=$NV/04_qc/compare_v2.4; mkdir -p $O; cd $O
V24=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/v2.4-20251210/01_final_results/01_assemblies/Manouria_impressa_v2.genome.fasta
minimap2 -x asm5 -t 48 --secondary=no $NV/06_release/Mimp_v3.fa $V24 > v2.4_to_v3.paf 2> minimap2.log
python3 - <<'PY'
import collections
hits=collections.defaultdict(collections.Counter); L={}; cov=collections.defaultdict(int)
for l in open("v2.4_to_v3.paf"):
    p=l.split("\t"); L[p[0]]=int(p[1])
    if int(p[10])<50000 or int(p[11])<30: continue
    hits[p[0]][p[5]]+=int(p[3])-int(p[2]); cov[p[0]]+=int(p[3])-int(p[2])
with open("v2.4_to_v3_correspondence.tsv","w") as f:
    f.write("v2.4_sequence\tlength\tv3_best_sequence\tbp_on_best\tpct_on_best\tsecond_v3_sequence\tpct_on_second\ttotal_aligned_pct\n")
    for s in sorted(L,key=lambda x:-L[x]):
        mc=hits[s].most_common(2); b=mc[0] if mc else ("-",0); c=mc[1] if len(mc)>1 else ("-",0)
        f.write(f"{s}\t{L[s]}\t{b[0]}\t{b[1]}\t{100*b[1]/L[s]:.1f}\t{c[0]}\t{100*c[1]/L[s]:.1f}\t{100*cov[s]/L[s]:.1f}\n")
split=[s for s in L if len(hits[s])>1 and hits[s].most_common(2)[1][1]/L[s]>0.05]
print("v2.4 sequences:",len(L),"; with >5% on a second v3 sequence (split/rearranged in v3):",len(split)); print("\n".join(split[:40]))
PY
head -n 40 v2.4_to_v3_correspondence.tsv; cp v2.4_to_v3_correspondence.tsv $NV/06_release/; echo "== DONE $(date)"
