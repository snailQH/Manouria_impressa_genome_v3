#!/bin/bash
#SBATCH --job-name=hificov
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/04c_hificov_%A.out
# Step 4c: HiFi read alignment to the scaffolds (minimap2 map-hifi) -> coordinate-sorted BAM -> depth in 10-kb
# windows. Used for: zero/low-coverage region detection (assembly artefacts), coverage-based contamination
# screening, mitochondrial contig identification, and the Technical Validation coverage figure.
#   usage: sbatch 04c_hifi_coverage.sbatch <label> <scaffolds.fa>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
LABEL=$1; FA=$(readlink -f $2); O=$NV/04_qc/hificov_$LABEL; mkdir -p $O; cd $O
echo "== $LABEL $(date) node=$(hostname)"; minimap2 --version
minimap2 -ax map-hifi -t 40 --secondary=no $FA $NV/00_data_qc/IT1.hifi.q20.filt.fastq.gz 2> minimap2.log \
  | samtools sort -@ 8 -m 3G -o hifi.bam - && samtools index -@ 8 hifi.bam
samtools flagstat -@ 8 hifi.bam > flagstat.txt; head -n 6 flagstat.txt
# mean depth per 10-kb window (MAPQ>=1 primary alignments)
samtools faidx $FA; awk '{print $1"\t"$2}' $FA.fai > sizes.txt
bedtools makewindows -g sizes.txt -w 10000 > win10k.bed
samtools view -@ 8 -b -q 1 -F 0x904 hifi.bam | bedtools coverage -a win10k.bed -b - -mean -sorted -g sizes.txt > depth10k.bed 2>/dev/null \
  || bedtools coverage -a win10k.bed -b <(samtools view -@ 8 -b -q 1 -F 0x904 hifi.bam) -mean > depth10k.bed
python3 - <<PY
import statistics, collections
d=[l.split() for l in open("depth10k.bed")]; vals=[float(x[3]) for x in d]
med=statistics.median(vals); print(f"windows {len(vals):,}; median depth {med:.1f}x; mean {statistics.mean(vals):.1f}x")
zero=[x for x in d if float(x[3])<0.05*med]; print(f"windows with depth <5% of median: {len(zero):,} ({sum(1 for x in zero if x[0].startswith('scaffold_') and int(x[0].split('_')[1])<=26):,} on chromosome-scale scaffolds)")
hi=[x for x in d if float(x[3])>2.5*med]; print(f"windows with depth >2.5x median (collapsed repeats / organelle): {len(hi):,}")
# merge zero windows into runs
runs=[]; 
for x in zero:
    if runs and runs[-1][0]==x[0] and int(x[1])==runs[-1][2]: runs[-1][2]=int(x[2])
    else: runs.append([x[0],int(x[1]),int(x[2])])
runs.sort(key=lambda r:-(r[2]-r[1]))
print("largest zero-coverage runs:"); [print(f"  {r[0]}:{r[1]}-{r[2]} ({(r[2]-r[1])/1e3:.0f} kb)") for r in runs[:15]]
open("zero_coverage_runs.bed","w").write("".join(f"{r[0]}\t{r[1]}\t{r[2]}\n" for r in runs))
PY
echo "== DONE $(date)"
