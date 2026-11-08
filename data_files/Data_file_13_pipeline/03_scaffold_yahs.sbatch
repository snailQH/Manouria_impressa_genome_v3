#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03_yahs_%A.out
# Step 3a: YaHS scaffolding from the name-sorted Hi-C BAM made by 03_hic_map.sbatch.
#   usage: sbatch 03_scaffold_yahs.sbatch <label> <contigs.fa> <hicmap label: primary|hap1|hap2> [extra yahs args]
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; CTG=$(readlink -f $2); MAPLABEL=$3; shift 3; EXTRA="$*"
O=$NV/03_scaffold/$LABEL; mkdir -p $O; cd $O
BAM0=$NV/03_scaffold/hicmap_$MAPLABEL/hic.namesorted.bam; [ -s $BAM0 ] || { echo "missing $BAM0"; exit 1; }
# chromap names mates NAME/1 and NAME/2; YaHS pairs consecutive records by identical QNAME, so strip the suffix.
BAM=$NV/03_scaffold/hicmap_$MAPLABEL/hic.namesorted.nosuffix.bam
[ -s $BAM ] || samtools view -h -@ 8 $BAM0 | awk 'BEGIN{OFS="\t"} /^@/{print; next} {sub(/\/[12]$/,"",$1); print}' | samtools view -@ 8 -b -o $BAM -
echo "== $LABEL $(date) yahs $(yahs --version) bam=$BAM extra='$EXTRA'"
yahs -o yahs -e GATC $EXTRA $CTG $BAM > yahs.log 2>&1
grep -E "read pairs" yahs.log | tail -n 1
assembly-stats -t yahs_scaffolds_final.fa | tee yahs_scaffolds_final.stats
python3 $NV/scripts/asm_stats.py yahs_scaffolds_final.fa | tee stats_full.tsv
echo "joins (U rows in AGP): $(awk '$5=="U"' yahs_scaffolds_final.agp | wc -l)"
# Juicebox files
juicer pre -a -o out_JBAT yahs.bin yahs_scaffolds_final.agp $CTG.fai > out_JBAT.log 2>&1
JT=/data1/reference/tools/juicer_tools.1.22.01.jar
java -Xmx100G -jar $JT pre out_JBAT.txt out_JBAT.hic.part <(grep PRE_C_SIZE out_JBAT.log | awk '{print $2" "$3}') > juicer_tools.log 2>&1 && mv out_JBAT.hic.part out_JBAT.hic
echo "== DONE $(date)"
