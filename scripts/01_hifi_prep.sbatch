#!/bin/bash
#SBATCH --job-name=hifi_prep
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/01_hifi_prep_%j.out
# Step 1: HiFi read preparation and k-mer survey.
#   - extract reads with CCS predicted accuracy rq >= 0.99 (Q20, the HiFi definition) from the Revio BAM
#   - screen for residual PacBio adapter (HiFiAdapterFilt-equivalent: BLAST of the SMRTbell adapter against reads is
#     replaced by an exact/near-exact search with seqkit; reads with a hit are removed and counted)
#   - meryl k=21 database and GenomeScope2 estimate of size/heterozygosity/repeats
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; O=$G/new-version/00_data_qc; cd $O
conda activate /data1/opt/conda/envs/assembly
echo "== versions"; samtools --version | head -1; seqkit version
echo "== extract rq>=0.99 $(date)"
samtools view -@ 16 -e '[rq]>=0.99' -b $G/01.pacbio_data/IT1.hifi_reads.bam | samtools fastq -@ 8 -T rq - | pigz -p 8 > IT1.hifi.q20.fastq.gz
echo "== adapter screen $(date)"
# PacBio SMRTbell adapter (v3) and its reverse complement; HiFiAdapterFilt uses BLAST at >=97% id over >=44 bp.
printf ">smrtbell_adapter\nATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT\n" > smrtbell.fa
seqkit grep -s -i -m 2 -p ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT -j 16 IT1.hifi.q20.fastq.gz | seqkit seq -n -i > adapter_hits.ids
echo "reads with adapter hit: $(wc -l < adapter_hits.ids)"
if [ -s adapter_hits.ids ]; then seqkit grep -v -f adapter_hits.ids -j 16 IT1.hifi.q20.fastq.gz -o IT1.hifi.q20.filt.fastq.gz; else cp IT1.hifi.q20.fastq.gz IT1.hifi.q20.filt.fastq.gz; fi
seqkit stats -a -j 16 -T IT1.hifi.q20.fastq.gz IT1.hifi.q20.filt.fastq.gz > hifi_stats.tsv; cat hifi_stats.tsv
echo "== meryl + genomescope2 $(date)"
conda activate /data1/opt/conda/envs/genomeqc
meryl k=21 count threads=32 memory=100 output IT1.k21.meryl IT1.hifi.q20.filt.fastq.gz
meryl histogram IT1.k21.meryl > IT1.k21.histo
genomescope2 -i IT1.k21.histo -o genomescope2_k21 -k 21 -p 2 > genomescope2_k21.log 2>&1 || true
cat genomescope2_k21/summary.txt 2>/dev/null
echo "== DONE $(date)"
