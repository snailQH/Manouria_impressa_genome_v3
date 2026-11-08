#!/bin/bash
#SBATCH --job-name=mask
#SBATCH --partition=main
#SBATCH --cpus-per-task=60
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/05b_mask_%A.out
# Step 12b: soft-mask the RELEASE assembly with the RepeatModeler library built in step 12 (05_repeats.sbatch)
# plus Dfam Tetrapoda; then classify the library's composition.  usage: sbatch 05b_mask_curated.sbatch <label> <genome.fa> <families.fa>
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
LABEL=$1; FA=$(readlink -f $2); LIB=$(readlink -f $3); O=$NV/05_annotation/$LABEL/01_repeats; mkdir -p $O; cd $O
SIF=/data1/reference/containers/tetools.sif; RUN="singularity exec -B /data2,/data1 $SIF"; export LC_ALL=C
echo "== $LABEL $(date) node=$(hostname)"; $RUN RepeatMasker -v; ln -sf $FA genome.fa; cp $LIB genome-families.fa
mkdir -p rm_custom rm_dfam
$RUN RepeatMasker -pa 15 -xsmall -gff -lib genome-families.fa -dir rm_custom genome.fa > rm_custom.log 2>&1
$RUN RepeatMasker -pa 15 -xsmall -gff -species tetrapoda -dir rm_dfam rm_custom/genome.fa.masked > rm_dfam.log 2>&1
# RepeatMasker uppercases its input: the Dfam-pass output carries only the Dfam masking. Union both passes:
grep -v '^#' rm_dfam/genome.fa.masked.out.gff | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' > rm_dfam/dfam_pass.bed
bedtools maskfasta -soft -fi rm_custom/genome.fa.masked -bed rm_dfam/dfam_pass.bed -fo genome.softmasked.fa
cat rm_custom/genome.fa.out.gff rm_dfam/genome.fa.masked.out.gff | grep -v "^#" | sort -k1,1 -k4,4n > repeats.gff
python3 - <<PY
L=0;m=0
for l in open("genome.softmasked.fa"):
    if l.startswith(">"): continue
    s=l.strip(); L+=len(s); m+=sum(1 for c in s if c.islower())
print(f"soft-masked {m:,} / {L:,} bp = {100*m/L:.2f}%")
PY
echo "== custom-library table"; cat rm_custom/genome.fa.tbl; echo "== Dfam pass table"; cat rm_dfam/genome.fa.masked.tbl
echo "== DONE $(date)"
