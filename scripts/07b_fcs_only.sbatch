#!/bin/bash
#SBATCH --job-name=fcs
#SBATCH --partition=main
#SBATCH --cpus-per-task=24
#SBATCH --mem=450G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/07b_fcs_%A.out
# Step 11b: FCS-adaptor + FCS-GX only (no MitoHiFi) — for the haplotype assemblies.  usage: sbatch 07b_fcs_only.sbatch <label> <fasta> <taxid>
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; C=/data1/reference/containers
LABEL=$1; FA=$(readlink -f $2); TAXID=$3; O=$NV/07_contamination/$LABEL; mkdir -p $O/adaptor $O/gx; cd $O
echo "== $LABEL $(date) node=$(hostname)"
singularity exec -B /data2,/data1 $C/fcs-adaptor.sif /app/fcs/bin/av_screen_x -o $O/adaptor --euk $FA > adaptor.log 2>&1; echo "adaptor: $(grep -vc "^#" $O/adaptor/fcs_adaptor_report.txt) flagged sequences"; cat $O/adaptor/fcs_adaptor_report.txt | head -n 20
GXDB=/data1/reference/fcs/gxdb
if [ -f $GXDB/.complete ]; then
  export FCS_DEFAULT_IMAGE=$C/fcs-gx.sif
  python3 /data1/reference/fcs/fcs.py --image $C/fcs-gx.sif screen genome --fasta $FA --out-dir $O/gx --gx-db $GXDB --tax-id $TAXID > gx.log 2>&1
  echo "gx: $(grep -vc '^#' $O/gx/*.fcs_gx_report.txt 2>/dev/null) flagged sequences"; head -n 20 $O/gx/*.fcs_gx_report.txt 2>/dev/null
else echo "gx: database not complete, skipped (rerun when /data1/reference/fcs/gxdb is complete)"; fi
echo "== DONE $(date)"
