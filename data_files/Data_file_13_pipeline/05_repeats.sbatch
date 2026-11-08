#!/bin/bash
#SBATCH --job-name=repeats
#SBATCH --partition=main
#SBATCH --cpus-per-task=60
#SBATCH --mem=200G
#SBATCH --time=168:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/05_repeats_%A.out
# Step 5: de novo repeat library and soft-masking (Dfam TETools container: RepeatModeler2 + RepeatMasker).
#   usage: sbatch 05_repeats.sbatch <label> <genome.fa>
#   RepeatModeler2 -LTRStruct builds the species library; RepeatMasker then masks with that library combined with
#   the Dfam curated Tetrapoda consensus set (-lib custom + -species via a second pass merged by union).
set -euo pipefail
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; FA=$(readlink -f $2); O=$NV/05_annotation/$LABEL/01_repeats; mkdir -p $O; cd $O
SIF=/data1/reference/containers/tetools.sif; RUN="singularity exec -B /data2,/data1 $SIF"
export LC_ALL=C
echo "== $LABEL $(date) node=$(hostname)"; $RUN RepeatModeler -version; $RUN RepeatMasker -v
ln -sf $FA genome.fa
[ -s genome.fa.nhr ] || $RUN BuildDatabase -name genome genome.fa > builddb.log 2>&1
if [ ! -s genome-families.fa ]; then $RUN RepeatModeler -database genome -threads 60 -LTRStruct > repeatmodeler.log 2>&1; fi
echo "families: $(grep -c '>' genome-families.fa)"
# pass 1: species-specific library; pass 2: Dfam Tetrapoda on the pass-1-masked genome; union via ProcessRepeats outputs
mkdir -p rm_custom rm_dfam
$RUN RepeatMasker -pa 15 -xsmall -gff -lib genome-families.fa -dir rm_custom genome.fa > rm_custom.log 2>&1
$RUN RepeatMasker -pa 15 -xsmall -gff -species tetrapoda -dir rm_dfam rm_custom/genome.fa.masked > rm_dfam.log 2>&1
cp rm_dfam/genome.fa.masked.masked genome.softmasked.fa
cat rm_custom/genome.fa.out.gff rm_dfam/genome.fa.masked.out.gff | grep -v "^#" | sort -k1,1 -k4,4n > repeats.gff
echo "== masked fraction"; python3 - <<PY
import re
L=0;m=0
for l in open("genome.softmasked.fa"):
    if l.startswith(">"): continue
    s=l.strip(); L+=len(s); m+=sum(1 for c in s if c.islower())
print(f"soft-masked {m:,} / {L:,} bp = {100*m/L:.2f}%")
PY
cat rm_custom/genome.fa.tbl | head -40
echo "== DONE $(date)"
