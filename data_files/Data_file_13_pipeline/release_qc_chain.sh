#!/bin/bash
# Waits for the release FASTA, then submits everything that depends only on it, across idle nodes:
# QC (compleasm/BUSCO), Merqury, contact map, synteny, HiFi coverage, FCS/mito, haplotype compleasm.
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; L=$1; FA=$NV/06_release/$L.fa
until [ -s $FA.fai ]; do sleep 60; done; sleep 30
S=$NV/scripts
Q=$(sbatch --parsable $S/04_qc.sbatch $L $FA)
M=$(sbatch --parsable $S/04b_merqury.sbatch $L $FA)
C=$(sbatch --parsable $S/03c_contactmap.sbatch $L $FA)
Y=$(sbatch --parsable $S/03d_synteny.sbatch $L $FA)
H=$(sbatch --parsable $S/04c_hifi_coverage.sbatch $L $FA)
F=$(sbatch --parsable $S/07_fcs_mito.sbatch $L $FA 74924)
A1=$(sbatch --parsable $S/04_qc.sbatch hap1_scaffolds $NV/03_scaffold/hap1_haphic/haphic_scaffolds.fa)
A2=$(sbatch --parsable $S/04_qc.sbatch hap2_scaffolds $NV/03_scaffold/hap2_haphic/haphic_scaffolds.fa)
echo "RELEASE-QC $L qc=$Q merqury=$M contactmap=$C synteny=$Y hificov=$H fcs_mito=$F hap1qc=$A1 hap2qc=$A2 $(date)" | tee -a $NV/logs/jobs.log
