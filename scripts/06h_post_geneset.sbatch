#!/bin/bash
#SBATCH --job-name=post_geneset
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06h_post_geneset_%A.out
# Step 14g (after 06g): GH18 chitinase loci vs the FINAL gene set (release IDs); per-chromosome gene counts and density for the
#   gene-summary figure; tRNA/rRNA in release GFF3 (ncRNA file with release names); manuscript-input copies.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; A=$NV/05_annotation/Mimp_v3; F=$A/06_final_geneset; R=$NV/06_release; O=$A/07_post_geneset; mkdir -p $O; cd $O
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/busco
diamond makedb --in $F/Mimp_v3.proteins.faa -d final_aa > makedb.log 2>&1 || { echo 'makedb failed'; cat makedb.log; exit 1; }
diamond blastp -q $NV/00_data_qc/gh18_queries.faa -d final_aa -o gh18_vs_final.tsv -f 6 qseqid sseqid pident length evalue bitscore --max-target-seqs 5 --quiet
echo "== GH18 queries vs final proteins (best hit per query)"; sort -k1,1 -k6,6nr gh18_vs_final.tsv | awk '!s[$1]++' | tee gh18_best_hits.tsv
echo "== final genes in the chr4 GH18 region (139.0–139.3 Mb)"; awk -F'\t' '$1=="chr4" && $3=="gene" && $4>139000000 && $5<139300000' $F/Mimp_v3.genes.gff3 | cut -f1,4,5,7,9 | tee gh18_region_genes.tsv
conda activate /data1/opt/conda/envs/assembly
awk -F'\t' '$3=="gene"{n[$1]++} END{for(k in n) print k"\t"n[k]}' $F/Mimp_v3.genes.gff3 | sort -k1,1V > genes_per_sequence.tsv
join -t $'\t' <(sort -k1,1 genes_per_sequence.tsv) <(cut -f1,2 $R/Mimp_v3.fa.fai | sort -k1,1) | awk -F'\t' 'BEGIN{OFS="\t"; print "sequence","genes","length_bp","genes_per_Mb"} {print $1,$2,$3,1e6*$2/$3}' | sort -k1,1V > gene_density.tsv
# ncRNA GFF with release names (tRNA high-confidence + rRNA) for deposit
N=$A/04_ncrna; awk -F'\t' 'NR>3 && $NF ~ /high confidence set/' $N/hc/trna_hc.out | awk 'NR>0{s=$3; e=$4; st="+"; if(s>e){t=s; s=e; e=t; st="-"}; printf "%s\ttRNAscan-SE\ttRNA\t%d\t%d\t%s\t%s\t.\tID=tRNA_%s_%s;product=tRNA-%s;anticodon=%s\n",$1,s,e,$9,st,$1,$2,$5,$6}' > Mimp_v3.tRNA_hc.gff3
grep -v "^#" $N/rrna.gff > Mimp_v3.rRNA.gff3 || true
echo "tRNA HC features: $(grep -c "" Mimp_v3.tRNA_hc.gff3); rRNA features: $(grep -c "" Mimp_v3.rRNA.gff3)"
cp Mimp_v3.tRNA_hc.gff3 Mimp_v3.rRNA.gff3 $R/; cp gene_density.tsv gh18_best_hits.tsv gh18_region_genes.tsv $NV/09_manuscript_inputs/
echo "== DONE $(date)"
