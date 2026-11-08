#!/bin/bash
# Writes 06_release/software_versions.tsv: every tool used by the v3 pipeline with its exact version (conda envs + containers). Table 1 source.
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; source /data1/opt/conda/etc/profile.d/conda.sh
OUT=$NV/06_release/software_versions.tsv; echo -e "tool\tversion\tenvironment\tstep" > $OUT
v(){ conda list -p /data1/opt/conda/envs/$1 2>/dev/null | awk -v t="$2" -v e="$1" -v s="$3" 'BEGIN{OFS="\t"} $1==t{print t,$2,e,s; f=1} END{if(!f) print t,"NOT FOUND",e,s}' >> $OUT; }
v assembly samtools "1,3"; v assembly seqkit "1,15"; v genomeqc meryl "1,4"; v genomeqc genomescope2 "1"; v assembly hifiasm "2"; v assembly gfatools "2"
v assembly chromap "3"; v assembly yahs "3"; v assembly minimap2 "3e,4c,15"; v assembly pretextmap "3c"; v assembly pretextsnapshot "3c"; v assembly bedtools "3f"
v busco busco "4,14"; v genomeqc merqury "4"; v assembly assembly-stats "4"; v busco diamond "14"; v eggnog eggnog-mapper "14"
echo -e "compleasm\t$(conda run -p /data1/opt/conda/envs/compleasm compleasm --version 2>/dev/null | grep -oE "[0-9.]+" | head -1) (bioconda package 0.2.6)\tcompleasm\t4" >> $OUT
echo -e "HapHiC\tcommit 74f1cb2\thaphic\t3b" >> $OUT; echo -e "juicer_tools\t1.22.01\tjar\t3" >> $OUT
C=/data1/reference/containers
echo -e "BRAKER3\t$(singularity exec $C/braker3.sif braker.pl --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)\tbraker3.sif\t13" >> $OUT
echo -e "RepeatModeler\t$(singularity exec $C/tetools.sif RepeatModeler -version 2>/dev/null | grep -oE '[0-9.]+' | head -1)\ttetools.sif\t5" >> $OUT
echo -e "RepeatMasker\t$(singularity exec $C/tetools.sif RepeatMasker -v 2>/dev/null | grep -oE '[0-9.]+[-a-z0-9]*' | head -1)\ttetools.sif\t12" >> $OUT
echo -e "FCS-adaptor\t$(singularity inspect $C/fcs-adaptor.sif 2>/dev/null | grep -oE 'fcs-adaptor:[0-9.]+' | cut -d: -f2 | head -1)\tfcs-adaptor.sif\t11" >> $OUT
echo -e "FCS-GX\t$(singularity exec $C/fcs-gx.sif gx --help 2>&1 | grep -oE 'git:v[0-9.]+' | cut -d: -f2 | head -1)\tfcs-gx.sif\t11" >> $OUT
echo -e "MitoHiFi\t$(singularity exec $C/mitohifi.sif mitohifi.py --version 2>/dev/null | head -1 || echo see mitohifi.log)\tmitohifi.sif\t11" >> $OUT
cat $OUT
