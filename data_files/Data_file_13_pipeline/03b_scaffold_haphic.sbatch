#!/bin/bash
#SBATCH --job-name=haphic
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03b_haphic_%A.out
# Step 3b: Hi-C scaffolding with HapHiC (prior peer review).  usage: sbatch 03b_scaffold_haphic.sbatch <label> <contigs.fa> <n_chromosomes>
#   Hi-C mapped with chromap (same as YaHS step) -> BAM (filtered as HapHiC recommends) -> haphic pipeline
#   -> scaffolds FASTA + AGP + Juicebox files.  n_chromosomes: 26 (Testudinidae 2n = 52).
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; CTG=$(readlink -f $2); NCHR=$3; MAPLABEL=$4; O=$NV/03_scaffold/$LABEL; mkdir -p $O; cd $O
HAPHIC=/data1/reference/tools/HapHiC; echo "== $LABEL $(date) HapHiC $(cd $HAPHIC && git log -1 --format=%h)"
samtools faidx $CTG
BAM=$NV/03_scaffold/hicmap_$MAPLABEL/hic.filtered.bam
ln -sf $BAM hic.filtered.bam; echo "bam=$BAM"
conda activate /data1/opt/conda/envs/haphic
$HAPHIC/haphic pipeline $CTG hic.filtered.bam $NCHR --threads 48 --processes 16 --RE GATC > haphic.log 2>&1 || { tail -n 20 haphic.log; exit 1; }
conda activate /data1/opt/conda/envs/assembly
cp 04.build/scaffolds.fa haphic_scaffolds.fa; cp 04.build/scaffolds.agp haphic_scaffolds.agp 2>/dev/null || true
assembly-stats -t haphic_scaffolds.fa | tee haphic_scaffolds.stats
# Juicebox files (HapHiC's build step writes a juicebox.sh with the juicer_tools command)
ls 04.build/ | head; cat 04.build/juicebox.sh 2>/dev/null | head -5
echo "== DONE $(date)"
