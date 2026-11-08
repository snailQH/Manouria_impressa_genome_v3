#!/bin/bash
#SBATCH --job-name=synteny_hic
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/16_synteny_hic_%A.out
# Step 15d: synteny blocks vs G. evgoodei, candidate structural differences, Hi-C support for each breakpoint in Mimp_v3.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; O=$NV/04_qc/synteny_rearrangements; mkdir -p $O; cd $O
python3 $NV/scripts/16_synteny_rearrangements.py $NV/03_scaffold/synteny_Mimp_v3/aln.paf $NV/03_scaffold/contactmap_Mimp_v3/hic.bam $NV/06_release/Mimp_v3.fa.fai $NV/06_release/Mimp_v3.names.tsv $NV/03_scaffold/curated/v3_curated.agp Mimp_v3_vs_GopEvg 1000000
cat Mimp_v3_vs_GopEvg_composition.tsv; echo; cat Mimp_v3_vs_GopEvg_breakpoints.tsv; cp Mimp_v3_vs_GopEvg_*.tsv $NV/09_manuscript_inputs/; echo "== DONE $(date)"
