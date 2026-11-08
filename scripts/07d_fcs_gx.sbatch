#!/bin/bash
#SBATCH --job-name=fcs_gx
#SBATCH --partition=main
#SBATCH --cpus-per-task=24
#SBATCH --mem=450G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/07d_fcs_gx_%A.out
# Step 11d: NCBI FCS-GX 0.5.5 only (gxdb complete flag required).  usage: sbatch 07d_fcs_gx.sbatch <label> <fasta> <taxid>
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; C=/data1/reference/containers; GXDB=/data1/reference/fcs/gxdb
LABEL=$1; FA=$(readlink -f $2); TAXID=$3; O=$NV/07_contamination/$LABEL/gx; mkdir -p $O; cd $O
[ -f $GXDB/.complete ] || { echo "gxdb not complete"; exit 1; }
echo "== $LABEL $(date) node=$(hostname)"
export FCS_DEFAULT_IMAGE=$C/fcs-gx.sif
python3 /data1/reference/fcs/fcs.py --image $C/fcs-gx.sif screen genome --fasta $FA --out-dir $O --gx-db $GXDB --tax-id $TAXID > gx.log 2>&1
echo "gx: $(grep -vc '^#' $O/*.fcs_gx_report.txt) flagged sequences"; cat $O/*.fcs_gx_report.txt | head -n 30; tail -n 3 gx.log
echo "== DONE $(date)"
