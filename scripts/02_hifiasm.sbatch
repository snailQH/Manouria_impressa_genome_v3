#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=main
#SBATCH --cpus-per-task=60
#SBATCH --mem=450G
#SBATCH --time=96:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/02_hifiasm_%j.out
# Step 2: contig assembly with hifiasm in Hi-C-integrated mode.
#   Produces a primary assembly and two phased haplotype assemblies (hap1/hap2) from HiFi + Hi-C,
#   as prior peer review requested. Default parameters otherwise; -l2 purge level is the hifiasm default in Hi-C mode.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; O=$G/new-version/01_hifiasm; cd $O
echo "== hifiasm $(hifiasm --version) $(date) node=$(hostname)"
hifiasm -o IT1 -t 60 --h1 $G/02.hi-c_data/IT-1_R1.fq.gz --h2 $G/02.hi-c_data/IT-1_R2.fq.gz $G/new-version/00_data_qc/IT1.hifi.q20.filt.fastq.gz 2> IT1.hifiasm.log
for g in IT1.hic.p_ctg IT1.hic.hap1.p_ctg IT1.hic.hap2.p_ctg; do
  gfatools gfa2fa $g.gfa > $g.fa; echo "== $g"; assembly-stats -t $g.fa
done
echo "== DONE $(date)"
