#!/bin/bash
#SBATCH --job-name=mito_verify
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/07e_mito_verify_%A.out
# Step 11e: verify the mitogenome structure with HiFi reads: depth along the molecule and the number of reads spanning the
#   tandem tRNA-Thr/tRNA-Pro repeat + control region (reference NC_011815.1 is 16,642 bp; MitoHiFi assembly 18,133 bp).
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; M=$NV/07_contamination/Mimp_v3/mito; cd $M
MT=final_mitogenome.fasta; L=$(awk '!/^>/{n+=length($0)} END{print n}' $MT)
# doubled reference so reads crossing the origin align contiguously
python3 - <<PY
s="".join(l.strip() for l in open("$MT") if not l.startswith(">")); open("mito_doubled.fa","w").write(">mito2x\n"+s+s+"\n")
PY
minimap2 -x map-hifi -t 16 -a mito_doubled.fa $NV/00_data_qc/IT1.hifi.q20.filt.fastq.gz 2>/dev/null | samtools view -b -F 2308 -q 20 - | samtools sort -@4 -o mito_reads.bam - && samtools index mito_reads.bam
echo "mito length $L; reads mapped (primary, MAPQ>=20): $(samtools view -c mito_reads.bam)"
samtools depth -a mito_reads.bam | awk -v L=$L '$2<=L{d[$2]=$3} END{n=asort(d); print "depth over molecule: median",d[int(n/2)],"min",d[1],"max",d[n]}'
# region of interest: tandem repeat + control region = from the last ND6/CYTB-side tRNA to the end (ref block ends at 15,547)
samtools view mito_reads.bam | awk -v L=$L '{pos=$4; cig=$6; len=0; while(match(cig,/[0-9]+[MDN=X]/)){len+=substr(cig,RSTART,RLENGTH-1); cig=substr(cig,RSTART+RLENGTH)} e=pos+len-1; if(pos<=15300 && e>=L+200) n++; if(pos<=15300 && e>=L) m++} END{print "reads spanning 15,300 -> past the origin (+200 bp):",n; print "reads spanning 15,300 -> origin:",m}'
samtools depth -a -r mito2x:15300-$L mito_reads.bam | awk '{s+=$3; n++} END{printf "mean depth 15,300-end: %.1f\n", s/n}'
echo "== DONE $(date)"
