#!/bin/bash
#SBATCH --job-name=verifycur
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03e_verify_%A.out
# Step 3e: verify curation — remap Hi-C to the curated scaffolds, re-score every join (join_qc.py), render
# Pretext map + per-scaffold snapshots.   usage: sbatch 03e_verify_curation.sbatch <label> <curated.fa> <curated.agp>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; FA=$(readlink -f $2); AGP=$(readlink -f $3); O=$NV/03_scaffold/contactmap_$LABEL; mkdir -p $O; cd $O
samtools faidx $FA; chromap -i -r $FA -o scaf.index
chromap --preset hic -r $FA -x scaf.index -q 1 --remove-pcr-duplicates -t 48 --SAM \
    -1 $G/02.hi-c_data/IT-1_R1.fq.gz -2 $G/02.hi-c_data/IT-1_R2.fq.gz -o hic.sam 2> chromap.log
samtools view -@ 16 -b -o hic.bam hic.sam && rm hic.sam
python3 $NV/scripts/join_qc.py hic.bam $AGP $FA.fai join_qc_$LABEL.tsv | tee join_qc_$LABEL.summary
samtools view -h hic.bam | PretextMap -o $LABEL.pretext --sortby length --mapq 1 > pretext.log 2>&1
PretextSnapshot -m $LABEL.pretext --sequences "=full" -r 2000 -o snapshots > /dev/null 2>&1 || true
PretextSnapshot -m $LABEL.pretext --sequences "=all" -r 1500 -o snapshots_each > /dev/null 2>&1 || true
echo "== DONE $(date)"
