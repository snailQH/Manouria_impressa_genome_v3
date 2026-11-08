#!/bin/bash
#SBATCH --job-name=braker3
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=168:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06_braker3_%A.out
# Step 6: gene prediction with BRAKER3 (GeneMark-ETP + AUGUSTUS) using RNA-seq (heart IT3, lung IT4) and the
# OrthoDB v12 Vertebrata protein set; input = soft-masked genome from step 5.
#   usage: sbatch 06_braker3.sbatch <label>
set -euo pipefail
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; NV=$G/new-version
LABEL=$1; A=$NV/05_annotation/$LABEL; O=$A/02_braker3; mkdir -p $O $A/augustus_config; cd $O
SIF=/data1/reference/containers/braker3.sif; RUN="singularity exec -B /data2,/data1 $SIF"
echo "== $LABEL $(date) node=$(hostname)"; $RUN braker.pl --version
rm -rf $A/augustus_config/species/Manouria_impressa_v3
# writable AUGUSTUS config copy (BRAKER needs to write species parameters)
[ -d $A/augustus_config/species ] || $RUN bash -c "cp -r \$AUGUSTUS_CONFIG_PATH/. $A/augustus_config/"
[ -s $NV/00_data_qc/orthodb_vertebrata.fa ] || zcat /data1/reference/orthodb/Vertebrata.fa.gz > $NV/00_data_qc/orthodb_vertebrata.fa
RNA=$G/03.RNA-Seq_data; mkdir -p $A/rnaseq
if [ ! -s $A/rnaseq/IT4.bam ]; then
  $RUN hisat2-build -p 48 $A/01_repeats/Mimp_v3.softmasked.fa $A/rnaseq/genome > $A/rnaseq/hisat2_build.log 2>&1
  for S in IT3 IT4; do
    $RUN bash -c "hisat2 -p 40 --dta -x $A/rnaseq/genome -1 $RNA/$S/${S}_R1.fq.gz -2 $RNA/$S/${S}_R2.fq.gz 2> $A/rnaseq/$S.hisat2.log | samtools sort -@ 8 -m 2G -o $A/rnaseq/$S.bam -"
    grep "overall alignment rate" $A/rnaseq/$S.hisat2.log
  done
fi
$RUN braker.pl --genome=$A/01_repeats/Mimp_v3.softmasked.fa --species=Manouria_impressa_v3 \
  --prot_seq=$NV/00_data_qc/orthodb_vertebrata.fa --bam=$A/rnaseq/IT3.bam,$A/rnaseq/IT4.bam \
  --AUGUSTUS_CONFIG_PATH=$A/augustus_config --threads=48 --workingdir=$O --gff3 > braker.log 2>&1
# (--busco_lineage is NOT used: BRAKER would run compleasm, which downloads its lineage from the internet; compute nodes have none.
#  Gene-set completeness is assessed afterwards by 06b_functional.sbatch with the offline BUSCO lineage.)
ls -la braker/braker.gff3 braker/braker.aa braker/braker.codingseq 2>/dev/null || ls -la braker.gff3 braker.aa braker.codingseq
G3=$(ls braker/braker.gff3 2>/dev/null || ls braker.gff3); echo "genes: $(grep -c $'\tgene\t' $G3)  mRNAs: $(grep -c $'\tmRNA\t' $G3)"
echo "== DONE $(date)"
