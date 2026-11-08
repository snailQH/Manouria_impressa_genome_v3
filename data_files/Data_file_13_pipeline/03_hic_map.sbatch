#!/bin/bash
#SBATCH --job-name=hicmap
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03_hicmap_%A.out
# Step 3.0: map Hi-C to a contig set ONCE, to disk (chromap --SAM to /dev/stdout aborts), and derive the
# inputs every scaffolder needs: name-sorted BAM (YaHS), filtered BAM (HapHiC), pairs stats.
#   usage: sbatch 03_hic_map.sbatch <label> <contigs.fa>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; CTG=$(readlink -f $2); O=$NV/03_scaffold/hicmap_$LABEL; mkdir -p $O; cd $O
echo "== $LABEL $(date) node=$(hostname)"; chromap --version 2>&1 | head -1
samtools faidx $CTG; ln -sf $CTG contigs.fa; ln -sf $CTG.fai contigs.fa.fai
[ -s contigs.index ] || chromap -i -r $CTG -o contigs.index
chromap --preset hic -r $CTG -x contigs.index -q 1 --remove-pcr-duplicates -t 48 --SAM \
    -1 $G/02.hi-c_data/IT-1_R1.fq.gz -2 $G/02.hi-c_data/IT-1_R2.fq.gz -o hic.sam 2> chromap.log
grep -E "Number of" chromap.log
samtools view -@ 16 -b -o hic.unsorted.bam hic.sam && rm hic.sam
samtools sort -@ 16 -n -m 3G -o hic.namesorted.bam hic.unsorted.bam          # YaHS input
samtools view -@ 16 -b -F 3340 -o hic.filtered.bam hic.unsorted.bam           # HapHiC input (their recommended filter)
rm hic.unsorted.bam
samtools flagstat -@ 8 hic.namesorted.bam > flagstat.txt; head -n 8 flagstat.txt
echo "== DONE $(date)"
