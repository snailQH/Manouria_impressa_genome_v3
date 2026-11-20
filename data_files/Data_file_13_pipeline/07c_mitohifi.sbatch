#!/bin/bash
#SBATCH --job-name=mitohifi
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=$PROJECT/01.Genome_Assembly/new-version/logs/07c_mitohifi_%A.out
# Step 11c: MitoHiFi 3.2.1 from HiFi reads with a pre-fetched reference (findMitoReference fails on compute nodes: no NCBI access).
#   Reference: Manouria impressa mitogenome NC_011815.1 (fetched on the head node with efetch -> 07_contamination/mito_ref/).
#   usage: sbatch 07c_mitohifi.sbatch <label>
set -euo pipefail
NV=$PROJECT/01.Genome_Assembly/new-version; C=$REFDIR/containers; L=$1
R=$NV/07_contamination/mito_ref; O=$NV/07_contamination/$L/mito; mkdir -p $O; cd $O
echo "== $L $(date) node=$(hostname)"; grep ">" $R/NC_011815.1.fasta
singularity exec -B <shared-storage>,<shared-storage> $C/mitohifi.sif mitohifi.py -r $NV/00_data_qc/IT1.hifi.q20.filt.fastq.gz -f $R/NC_011815.1.fasta -g $R/NC_011815.1.gb -t 16 -o 2 > mitohifi.log 2>&1 || true
ls; (grep ">" final_mitogenome.fasta || true); tail -n 5 mitohifi.log
[ -s final_mitogenome.fasta ] && { source $CONDA_ROOT/etc/profile.d/conda.sh; conda activate $CONDA_ROOT/envs/assembly; seqkit stats final_mitogenome.fasta; }
echo "== DONE $(date)"
