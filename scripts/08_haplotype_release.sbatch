#!/bin/bash
#SBATCH --job-name=hap_release
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/08_hap_release_%A.out
# Step 15b: prepare the two haplotype assemblies for deposit as alternate assemblies:
#   rename scaffolds by length (h1_scaffold_N / h2_scaffold_N), drop < 1 kb, write release FASTA + stats,
#   align each haplotype to the primary release (minimap2 asm5) -> hapN_to_primary.tsv (best chr per scaffold).
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; R=$NV/06_release; cd $R
for H in 1 2; do
  IN=$NV/03_scaffold/hap${H}_haphic/haphic_scaffolds.fa; OUT=Mimp_v3_hap$H
  seqkit sort -l -r $IN 2>/dev/null | seqkit seq -m 1000 2>/dev/null | seqkit replace -p '.*' -r "h${H}_scaffold_{nr}" 2>/dev/null > $OUT.fa
  seqkit fx2tab -n -l $IN | sort -k2,2nr | awk -v h=$H 'BEGIN{OFS="\t"; print "release_name","haphic_name","length"} $2>=1000{print "h"h"_scaffold_"NR,$1,$2}' > $OUT.names.tsv
  samtools faidx $OUT.fa; python3 $NV/scripts/asm_stats.py $OUT.fa > $OUT.stats.tsv; cat $OUT.stats.tsv | head -n 6
  minimap2 -x asm5 -t 48 --secondary=no Mimp_v3.fa $OUT.fa > ${OUT}_to_primary.paf 2> /dev/null
  python3 - $OUT <<PY
import sys, collections
o=sys.argv[1]; hits=collections.defaultdict(collections.Counter); L={}
for l in open(o+"_to_primary.paf"):
    p=l.split("\t"); L[p[0]]=int(p[1])
    if int(p[10])<20000 or int(p[11])<30: continue
    hits[p[0]][p[5]]+=int(p[3])-int(p[2])
with open(o+"_to_primary.tsv","w") as f:
    f.write("hap_scaffold\tlength\tbest_primary_sequence\taligned_bp\tpct_of_scaffold\n")
    for s in sorted(L,key=lambda x:-L[x]):
        b=hits[s].most_common(1)[0] if hits[s] else ("-",0); f.write(f"{s}\t{L[s]}\t{b[0]}\t{b[1]}\t{100*b[1]/L[s]:.1f}\n")
chrom=sum(1 for s in L if hits[s] and hits[s].most_common(1)[0][0].startswith("chr") and hits[s].most_common(1)[0][1]/L[s]>0.5)
print(o, ": scaffolds", len(L), "; assigned to a chromosome (>50% aligned):", chrom)
PY
done
md5sum Mimp_v3.fa Mimp_v3.agp Mimp_v3_hap1.fa Mimp_v3_hap2.fa > release_md5.txt; cat release_md5.txt
echo "== DONE $(date)"
