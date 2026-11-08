#!/bin/bash
#SBATCH --job-name=pretext_render
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/14_pretext_render_%A.out
# Step 15c: publication renders of the Hi-C maps (PretextSnapshot, 2160 px, colour map "Red 2") for the release and both haplotypes.
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; O=$NV/09_manuscript_inputs; mkdir -p $O
for L in Mimp_v3 Mimp_v3_hap1 Mimp_v3_hap2; do
  M=$NV/03_scaffold/contactmap_$L/$L.pretext; [ -s $M ] || { echo "no map for $L"; continue; }
  PretextSnapshot -m $M --sequences "=full" -r 2160 -c "Red 2" -o $NV/03_scaffold/contactmap_$L/snapshots_pub --prefix ${L}_ > /dev/null 2>&1 || PretextSnapshot -m $M --sequences "=full" -r 2160 -c "Red 2" -o $NV/03_scaffold/contactmap_$L/snapshots_pub > /dev/null 2>&1
  f=$(ls $NV/03_scaffold/contactmap_$L/snapshots_pub/*.png | head -1); cp "$f" $O/hic_map_${L}_pub.png; echo "$L → $O/hic_map_${L}_pub.png"
done
echo "== DONE $(date)"
