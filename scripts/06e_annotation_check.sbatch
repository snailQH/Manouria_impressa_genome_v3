#!/bin/bash
#SBATCH --job-name=annot_check
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06e_annot_check_%A.out
# Step 14d: (1) BUSCO-protein of the component gene sets (GeneMark-ETP, AUGUSTUS) and of a "braker + non-overlapping GeneMark-ETP" union,
#   to decide whether the conservative TSEBRA set (C 89.6 %) should be supplemented; (2) RNA-seq read support per gene; (3) chitinase
#   (GH18) locus check: known M. impressa GH18 proteins mapped to Mimp_v3 and searched against the gene set.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; A=$NV/05_annotation/Mimp_v3; B=$A/02_braker3; O=$A/05_annotation_check; mkdir -p $O; cd $O
G=$A/01_repeats/Mimp_v3.softmasked.fa; SIF=/data1/reference/containers/braker3.sif; RUN="singularity exec -B /data2,/data1 $SIF"
source /data1/opt/conda/etc/profile.d/conda.sh
# ---- component sets → proteins
$RUN /opt/ETP/tools/gffread -y genemark.faa -g $G $B/GeneMark-ETP/genemark.gtf 2> gffread_gm.log || true
$RUN /opt/ETP/tools/gffread -y augustus.faa -g $G $B/Augustus/augustus.hints.gtf 2> gffread_aug.log || true
# ---- union: braker genes + GeneMark-ETP genes whose span does not overlap any braker gene (same strand)
conda activate /data1/opt/conda/envs/assembly
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5,$9,".",$7}' $B/braker.gtf | sort -k1,1 -k2,2n > braker_genes.bed
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5,$9,".",$7}' $B/GeneMark-ETP/genemark.gtf | sort -k1,1 -k2,2n > genemark_genes.bed
bedtools intersect -v -s -a genemark_genes.bed -b braker_genes.bed > genemark_extra_genes.bed
cut -f4 genemark_extra_genes.bed | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort -u > genemark_extra_ids.txt
echo "braker genes: $(wc -l < braker_genes.bed); GeneMark-ETP genes: $(wc -l < genemark_genes.bed); non-overlapping GeneMark-ETP genes: $(wc -l < genemark_extra_ids.txt)"
grep -F -w -f genemark_extra_ids.txt $B/GeneMark-ETP/genemark.gtf | sed 's/gene_id "\([^"]*\)"/gene_id "gm_\1"/; s/transcript_id "\([^"]*\)"/transcript_id "gm_\1"/' > genemark_extra.gtf
cat $B/braker.gtf genemark_extra.gtf > braker_plus_etp.gtf
$RUN /opt/ETP/tools/gffread -y braker_plus_etp.faa -g $G braker_plus_etp.gtf 2> gffread_plus.log || true
# longest isoform per gene for BUSCO
for S in genemark augustus braker_plus_etp; do python3 - $S.faa $S.longest.faa <<'PY'
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
for S in genemark augustus braker_plus_etp; do
  busco -i $S.longest.faa -l vertebrata_odb10 -m proteins -c 32 -o busco_$S --download_path /data1/reference/busco_downloads/busco_downloads --offline > busco_$S.log 2>&1 || true
  echo "== $S: $(grep -h 'C:' busco_$S/short_summary*.txt)"
done
echo "== braker (TSEBRA): $(grep -h 'C:' $A/03_functional/busco_prot_vertebrata/short_summary*.txt 2>/dev/null)"
# ---- RNA-seq support: reads per gene span (braker set), heart + lung
conda activate /data1/opt/conda/envs/assembly
for S in IT3 IT4; do [ -s $A/rnaseq/$S.bam.bai ] || samtools index $A/rnaseq/$S.bam; done
bedtools multicov -bams $A/rnaseq/IT3.bam $A/rnaseq/IT4.bam -bed braker_genes.bed > braker_genes_rnaseq_counts.bed
awk '{t++; if($7+$8>=10) s++; if($7+$8>=1) a++} END{printf "RNA-seq support (braker genes): >=10 reads %d/%d (%.1f %%); >=1 read %d (%.1f %%)\n", s,t,100*s/t,a,100*a/t}' braker_genes_rnaseq_counts.bed
# ---- chitinase locus check
conda activate /data1/opt/conda/envs/annot
miniprot -t 32 --gff $NV/06_release/Mimp_v3.fa $NV/00_data_qc/gh18_queries.faa > gh18_miniprot.gff 2> miniprot.log
awk -F'\t' '$3=="mRNA"{match($9,/Target=[^ ;]+/); t=substr($9,RSTART+7,RLENGTH-7); match($9,/Identity=[0-9.]+/); id=substr($9,RSTART+9,RLENGTH-9); print t, $1, $4, $5, $7, "identity="id}' gh18_miniprot.gff | sort -k2,2 -k3,3n | awk '$6>="identity=0.80"' > gh18_loci_on_Mimp_v3.txt
echo "== GH18 query placements (identity >= 0.80):"; cat gh18_loci_on_Mimp_v3.txt | head -40
conda activate /data1/opt/conda/envs/busco
diamond makedb --in $B/braker.aa -d braker_aa > /dev/null 2>&1; diamond blastp -q $NV/00_data_qc/gh18_queries.faa -d braker_aa -o gh18_vs_braker.tsv -f 6 qseqid sseqid pident length evalue bitscore --max-target-seqs 3 --quiet
echo "== GH18 queries vs braker proteins (best hit):"; sort -k1,1 -k6,6nr gh18_vs_braker.tsv | awk '!s[$1]++' | head -30
diamond makedb --in braker_plus_etp.faa -d plus_aa > /dev/null 2>&1; diamond blastp -q $NV/00_data_qc/gh18_queries.faa -d plus_aa -o gh18_vs_plus.tsv -f 6 qseqid sseqid pident length evalue bitscore --max-target-seqs 3 --quiet
echo "== GH18 queries vs braker+ETP proteins (best hit):"; sort -k1,1 -k6,6nr gh18_vs_plus.tsv | awk '!s[$1]++' | head -30
echo "== DONE $(date)"
