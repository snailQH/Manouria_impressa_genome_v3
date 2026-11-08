#!/bin/bash
#SBATCH --job-name=hic_stats
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03f_hic_stats_%A.out
# Step 3f: Hi-C library and mapping statistics for Table 2 / Technical Validation (prior peer review: "numbers and percentages of valid pairs").
#   Sources: chromap.log (mapping), hic.filtered.bam (dedup, MAPQ>=1, proper pairs): cis/trans and cis>=20 kb fractions.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; cd $NV/03_scaffold/hicmap_primary
TOT=$(grep "^Number of reads:" chromap.log | grep -oE "[0-9]+"); MAP=$(grep "^Number of mapped reads:" chromap.log | grep -oE "[0-9]+")
UNI=$(grep "^Number of uniquely mapped reads:" chromap.log | grep -oE "[0-9]+"); OUT=$(grep "passed filters" chromap.log | grep -oE "[0-9]+$")
samtools view -@16 -f 64 hic.filtered.bam | awk -v tot=$TOT -v map=$MAP -v uni=$UNI -v out=$OUT 'BEGIN{OFS="\t"}
 {n++; if($7=="="){cis++; d=$9<0?-$9:$9; if(d>=20000)far++; if(d>=1000)mid++} else trans++}
 END{print "metric","value","percent_of_pairs";
  print "read_pairs_sequenced",tot/2,"100";
  printf "%s\t%d\t%.1f\n","reads_mapped",map,100*map/tot; printf "%s\t%d\t%.1f\n","reads_uniquely_mapped",uni,100*uni/tot;
  printf "%s\t%d\t%.1f\n","pairs_after_dedup_MAPQ1_proper (valid pairs used for scaffolding)",n,200*n/tot;
  printf "%s\t%d\t%.1f\n","valid_pairs_cis",cis,100*cis/n; printf "%s\t%d\t%.1f\n","valid_pairs_cis_ge1kb",mid,100*mid/n;
  printf "%s\t%d\t%.1f\n","valid_pairs_cis_ge20kb",far,100*far/n; printf "%s\t%d\t%.1f\n","valid_pairs_trans",trans,100*trans/n;
  print "note","cis/trans measured on primary contigs (460 contigs), so trans includes intra-chromosomal inter-contig pairs",""}' > hic_stats.tsv
cat hic_stats.tsv; cp hic_stats.tsv $NV/06_release/Mimp_v3.hic_stats.tsv; echo "== DONE $(date)"
