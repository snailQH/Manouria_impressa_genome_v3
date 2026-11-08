#!/bin/bash
# Waits for the FCS-GX database download (curl, head node) to finish and verifies every file size against the NCBI listing,
# writes gxdb/.complete, then submits 07d_fcs_gx for primary, hap1, hap2 (sequential: each needs ~450 GB RAM).
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; D=/data1/reference/fcs/gxdb; S=$NV/scripts
URL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/
while pgrep -f "curl -sL -o gxdb/all" > /dev/null; do sleep 300; done
ok=1; for f in $(curl -sL $URL | grep -oE 'href="all\.[A-Za-z0-9._]+"' | cut -d'"' -f2 | sort -u); do
  r=$(curl -sIL $URL$f | grep -i "^content-length" | tail -1 | tr -dc 0-9); l=$(stat -c %s $D/$f 2>/dev/null || echo 0)
  [ "$r" = "$l" ] || { echo "MISMATCH $f remote=$r local=$l"; ok=0; }; done
if [ $ok = 1 ]; then touch $D/.complete; echo "gxdb complete $(date)"
  A=$(sbatch --parsable $S/07d_fcs_gx.sbatch Mimp_v3 $NV/06_release/Mimp_v3.fa 74924)
  B=$(sbatch --parsable --dependency=afterany:$A $S/07d_fcs_gx.sbatch Mimp_v3_hap1 $NV/06_release/Mimp_v3_hap1.fa 74924)
  C=$(sbatch --parsable --dependency=afterany:$B $S/07d_fcs_gx.sbatch Mimp_v3_hap2 $NV/06_release/Mimp_v3_hap2.fa 74924)
  echo "FCS-GX primary=$A hap1=$B hap2=$C $(date)" | tee -a $NV/logs/jobs.log
else echo "gxdb INCOMPLETE — re-download the mismatched files, then rerun this script" | tee -a $NV/logs/jobs.log; fi
