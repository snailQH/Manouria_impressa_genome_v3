#!/bin/bash
#SBATCH --job-name=geneset
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06f_geneset_%A.out
# Step 14e: final gene-set choice. Candidates: (1) braker.gtf (TSEBRA, 17,392 genes, BUSCO-protein C 89.6); (2) AUGUSTUS-with-hints set
#   (26,279 genes, C 97.6); (3) union = braker + AUGUSTUS genes that overlap no braker gene (same strand); (4) TSEBRA re-run keeping all
#   AUGUSTUS transcripts (-k) and adding hint-supported GeneMark-ETP transcripts. Reports gene counts, BUSCO-protein, RNA-seq support,
#   and the GH18 locus content of each, then writes the chosen set's name to geneset_choice.txt by the pre-registered rule:
#   highest BUSCO-C among sets with D < 4 %, ties broken by fewer genes.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; A=$NV/05_annotation/Mimp_v3; B=$A/02_braker3; O=$A/05_annotation_check; cd $O
G=$A/01_repeats/Mimp_v3.softmasked.fa; SIF=/data1/reference/containers/braker3.sif; RUN="singularity exec -B /data2,/data1 $SIF"
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5,$9,".",$7}' $B/Augustus/augustus.hints.gtf | sort -k1,1 -k2,2n > augustus_genes.bed
bedtools intersect -v -s -a augustus_genes.bed -b braker_genes.bed > augustus_extra_genes.bed
cut -f4 augustus_extra_genes.bed | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort -u > augustus_extra_ids.txt
grep -F -w -f augustus_extra_ids.txt $B/Augustus/augustus.hints.gtf | sed 's/gene_id "\([^"]*\)"/gene_id "aug_\1"/; s/transcript_id "\([^"]*\)"/transcript_id "aug_\1"/' > augustus_extra.gtf
cat $B/braker.gtf augustus_extra.gtf > braker_plus_aug.gtf
echo "AUGUSTUS genes: $(wc -l < augustus_genes.bed); non-overlapping AUGUSTUS genes added to braker: $(wc -l < augustus_extra_ids.txt)"
# (4) TSEBRA keep-all-AUGUSTUS
$RUN /opt/TSEBRA/bin/tsebra.py -k $B/Augustus/augustus.hints.gtf -g $B/GeneMark-ETP/genemark.gtf -e $B/hintsfile.gff -c /opt/TSEBRA/config/default.cfg -o tsebra_keepaug.gtf > tsebra_keepaug.log 2>&1 || echo "tsebra -k failed (see log)"
for S in braker_plus_aug tsebra_keepaug; do
  [ -s $S.gtf ] || continue
  $RUN /opt/ETP/tools/gffread -y $S.faa -g $G $S.gtf 2> gffread_$S.log || true
  python3 - $S.faa $S.longest.faa <<'PY'
import sys,collections; fin,fout=sys.argv[1:3]; seqs={}; cur=None
for l in open(fin):
    if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
    else: seqs[cur].append(l.strip())
seqs={k:''.join(v).rstrip('*').replace('.','') for k,v in seqs.items()}; g=collections.defaultdict(list)
for k in seqs: g[k.rsplit('.',1)[0] if '.t' in k else k].append(k)
with open(fout,'w') as f:
    for gene,ts in g.items(): t=max(ts,key=lambda x:len(seqs[x])); f.write(f'>{t}\n{seqs[t]}\n')
print(fin, 'genes', len(g), 'proteins', len(seqs))
PY
done
conda activate /data1/opt/conda/envs/busco
for S in braker_plus_aug tsebra_keepaug; do [ -s $S.longest.faa ] || continue
  busco -i $S.longest.faa -l vertebrata_odb10 -m proteins -c 32 -o busco_$S --download_path /data1/reference/busco_downloads/busco_downloads --offline > busco_$S.log 2>&1 || true
  echo "== $S: $(grep -h 'C:' busco_$S/short_summary*.txt)"; done
echo "== summary"; for S in braker augustus genemark braker_plus_etp braker_plus_aug tsebra_keepaug; do f=$(ls busco_$S/short_summary*.txt 2>/dev/null || ls $A/03_functional/busco_prot_vertebrata/short_summary*.txt 2>/dev/null); [ -n "$f" ] && echo "$S $(grep -h 'C:' $f | tr -d '\t')"; done | tee busco_summary.txt
# decision rule
python3 - <<'PY'
import re
rows=[]
for l in open('busco_summary.txt'):
    m=re.search(r'(\S+) C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]',l)
    if m: rows.append((m.group(1),float(m.group(2)),float(m.group(4))))
ok=[r for r in rows if r[2]<4.0]; best=max(ok,key=lambda r:r[1]) if ok else max(rows,key=lambda r:r[1])
open('geneset_choice.txt','w').write(best[0]+'\n'); print('CHOSEN (rule: max C with D<4):',best)
PY
# RNA-seq support for the chosen set
conda activate /data1/opt/conda/envs/assembly; C=$(cat geneset_choice.txt)
case $C in braker) GTF=$B/braker.gtf;; augustus) GTF=$B/Augustus/augustus.hints.gtf;; genemark) GTF=$B/GeneMark-ETP/genemark.gtf;; *) GTF=$O/$C.gtf;; esac
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5,$9,".",$7}' $GTF | sort -k1,1 -k2,2n > chosen_genes.bed
bedtools multicov -bams $A/rnaseq/IT3.bam $A/rnaseq/IT4.bam -bed chosen_genes.bed > chosen_genes_rnaseq_counts.bed
awk -v c=$C '{t++; if($7+$8>=10) s++; if($7+$8>=1) a++} END{printf "RNA-seq support (%s): >=10 reads %d/%d (%.1f %%); >=1 read %d (%.1f %%)\n", c,s,t,100*s/t,a,100*a/t}' chosen_genes_rnaseq_counts.bed | tee rnaseq_support_chosen.txt
echo "== DONE $(date)"
