#!/bin/bash
#SBATCH --job-name=fcs_mito
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=450G
#SBATCH --time=48:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/07_fcs_mito_%A.out
# Step 11: contamination screen (NCBI FCS-adaptor + FCS-GX) and mitochondrial genome (MitoHiFi) for the release assembly.
#   usage: sbatch 07_fcs_mito.sbatch <label> <release.fa> <taxid>      (M. impressa taxid 74924)
set -uo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; C=/data1/reference/containers
LABEL=$1; FA=$(readlink -f $2); TAXID=$3; O=$NV/07_contamination/$LABEL; mkdir -p $O/adaptor $O/gx $O/mito; cd $O
echo "== $LABEL $(date) node=$(hostname)"
# --- FCS-adaptor (vector/adapter) ---
singularity exec -B /data2,/data1 $C/fcs-adaptor.sif /app/fcs/bin/av_screen_x -o $O/adaptor --euk $FA > adaptor.log 2>&1; echo "adaptor: $(grep -vc "^#" $O/adaptor/fcs_adaptor_report.txt) flagged sequences"; cat $O/adaptor/fcs_adaptor_report.txt | head -n 20
# --- FCS-GX (foreign organisms); needs the gxdb (~470 GB) on local disk ---
GXDB=/data1/reference/fcs/gxdb
if [ -f $GXDB/.complete ]; then
  export FCS_DEFAULT_IMAGE=$C/fcs-gx.sif
  python3 /data1/reference/fcs/fcs.py --image $C/fcs-gx.sif screen genome --fasta $FA --out-dir $O/gx --gx-db $GXDB --tax-id $TAXID > gx.log 2>&1
  echo "gx: $(grep -vc '^#' $O/gx/*.fcs_gx_report.txt 2>/dev/null) flagged sequences"; head -n 20 $O/gx/*.fcs_gx_report.txt 2>/dev/null
else echo "gx: database not present, skipped (rerun when /data1/reference/fcs/gxdb is complete)"; fi
# --- MitoHiFi from HiFi reads, reference = closest chelonian mitogenome fetched by findMitoReference ---
cd $O/mito
singularity exec -B /data2,/data1 $C/mitohifi.sif findMitoReference.py --species "Manouria impressa" --outfolder . --min_length 15000 > findref.log 2>&1 || true
REF=$(ls *.fasta 2>/dev/null | head -1); GB=$(ls *.gb 2>/dev/null | head -1)
if [ -n "$REF" ]; then
  singularity exec -B /data2,/data1 $C/mitohifi.sif mitohifi.py -r $NV/00_data_qc/IT1.hifi.q20.filt.fastq.gz -f $REF -g $GB -t 48 -o 2 > mitohifi.log 2>&1
  ls; grep -iE "length|circular" final_mitogenome.fasta 2>/dev/null | head -2
else echo "mito: no reference fetched"; fi
echo "== DONE $(date)"
