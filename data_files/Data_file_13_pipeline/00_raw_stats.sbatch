#!/bin/bash
#SBATCH --job-name=raw_stats
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/raw_stats_%j.out
set -uo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/chitin_meta
G=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly; O=$G/new-version/00_data_qc
echo "== HiFi $(date)"
samtools view -@ 8 $G/01.pacbio_data/IT1.hifi_reads.bam | awk -F"\t" "{L=length(\$10); n++; b+=L; rq=\"NA\"; for(i=12;i<=NF;i++) if(\$i ~ /^rq:f:/) rq=substr(\$i,6); mv=substr(\$1,1,index(\$1,\"/\")-1); print L\"\t\"rq\"\t\"mv}" > $O/hifi_read_lengths.tsv
python3 - <<PY
import collections
L=[];q=[];mv=collections.Counter();lowq=0
for l in open("$O/hifi_read_lengths.tsv"):
    a,b,c=l.rstrip().split("\t"); n=int(a); L.append(n); mv[c]+=n
    if b!="NA" and float(b)<0.99: lowq+=1
L.sort(reverse=True); t=sum(L); a=0
for x in L:
    a+=x
    if a>=t/2: n50=x; break
print(f"HiFi reads {len(L):,} bases {t:,} ({t/1e9:.2f} Gb) mean {t/len(L):.0f} N50 {n50} max {L[0]} coverage@2.43Gb {t/2.43e9:.1f}x reads_rq<0.99 {lowq:,}")
for k,v in sorted(mv.items()): print("  movie",k,f"{v/1e9:.2f} Gb")
PY
echo "== Hi-C $(date)"
for r in R1 R2; do echo -n "$r pairs: "; zcat $G/02.hi-c_data/IT-1_$r.fq.gz | awk "NR%4==2{n++; b+=length(\$0)} END{print n\" reads, \"b\" bases, mean \"b/n}"; done
echo "== RNA $(date)"
for s in IT3 IT4; do echo -n "$s: "; zcat $G/03.RNA-Seq_data/$s/${s}_R1.fq.gz | awk "NR%4==2{n++; b+=length(\$0)} END{print n\" pairs, \"2*b\" bases\"}"; done
echo "== DONE $(date)"
