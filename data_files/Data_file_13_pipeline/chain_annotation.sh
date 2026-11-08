#!/bin/bash
# Waits for the RepeatModeler library, then submits masking of the RELEASE assembly -> BRAKER3 with dependencies.
# Run once with nohup on cluster.   usage: chain_annotation.sh <release label> <release.fa>
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; LABEL=$1; FA=$2
LIB=$NV/05_annotation/v3_yahs/01_repeats/genome-families.fa
until [ -s $LIB ] && grep -q "families" $NV/05_annotation/v3_yahs/01_repeats/repeatmodeler.log && ! pgrep -f "RepeatModeler -database" > /dev/null; do sleep 300; done
echo "library ready $(date): $(grep -c '>' $LIB) families"
M=$(sbatch --parsable $NV/scripts/05b_mask_curated.sbatch $LABEL $FA $LIB)
B=$(sbatch --parsable --dependency=afterok:$M $NV/scripts/06_braker3.sbatch $LABEL)
echo "ANNOTATION $LABEL mask=$M braker3=$B" | tee -a $NV/logs/jobs.log
