#!/bin/bash
#SBATCH -J gwh_fix
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -o /data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/06_release/gwh_fix/gwh_fix.%j.log
set -euo pipefail
cd /data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/06_release; mkdir -p gwh_fix
python3 ../scripts/06i_gff3_gwh_fix.py Mimp_v3.fa Mimp_v3.genes.gff3 Mimp_v3.proteins.faa Mimp_v3.cds.fna gwh_fix/Mimp_v3.genes
echo DONE
