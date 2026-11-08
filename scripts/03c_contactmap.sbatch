#!/bin/bash
#SBATCH --job-name=contactmap
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03c_contactmap_%A.out
# Step 3c: Hi-C contact map of a scaffold set (Pretext), for inspection, curation and Figure 2C.
#   usage: sbatch 03c_contactmap.sbatch <label> <scaffolds.fa>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; FA=$(readlink -f $2); O=$NV/03_scaffold/contactmap_$LABEL; mkdir -p $O; cd $O
samtools faidx $FA
chromap -i -r $FA -o scaf.index
chromap --preset hic -r $FA -x scaf.index -q 1 --remove-pcr-duplicates -t 48 --SAM \
    -1 $G/02.hi-c_data/IT-1_R1.fq.gz -2 $G/02.hi-c_data/IT-1_R2.fq.gz -o hic.sam 2> chromap.log
samtools view -@ 16 -b -o hic.bam hic.sam && rm hic.sam
samtools view -h hic.bam | PretextMap -o $LABEL.pretext --sortby length --mapq 1 > pretext.log 2>&1
PretextSnapshot -m $LABEL.pretext --sequences "=full" -r 2000 -o snapshots > /dev/null 2>&1 || true
ls snapshots; echo "== DONE $(date)"
