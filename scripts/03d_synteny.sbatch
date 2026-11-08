#!/bin/bash
#SBATCH --job-name=synteny
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/03d_synteny_%A.out
# Step 3d: whole-genome alignment of the v3 scaffolds to Gopherus evgoodei rGopEvg1 (GCF_007399415.2, VGP,
# chromosome-level) with minimap2 -x asm20; per-scaffold best-hit chromosome assignment and a dot plot.
#   usage: sbatch 03d_synteny.sbatch <label> <scaffolds.fa>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version
LABEL=$1; FA=$(readlink -f $2); O=$NV/03_scaffold/synteny_$LABEL; mkdir -p $O; cd $O
REF=/data1/reference/genomes/GCF_007399415.2_rGopEvg1_v1.p_genomic.fna.gz
# keep only chromosome-assigned reference sequences, renamed to their chromosome label
[ -s ref_chr.fa ] || python3 - <<PY
import gzip,re
out=open("ref_chr.fa","w"); keep=False; n=0
for l in gzip.open("$REF","rt"):
    if l.startswith(">"):
        m=re.search(r"chromosome (\w+),", l); keep=bool(m) and "unlocalized" not in l
        if keep: out.write(">GopEvg_chr"+m.group(1)+" "+l[1:].split()[0]+"\n"); n+=1
    elif keep: out.write(l)
print("reference chromosomes kept:",n)
PY
grep -c ">" ref_chr.fa
minimap2 -x asm20 -t 48 --secondary=no ref_chr.fa $FA > aln.paf 2> minimap2.log
python3 - <<PY
import collections
# per query scaffold: aligned bp per reference chromosome (alignments >= 50 kb, mapq >= 30)
hits=collections.defaultdict(lambda: collections.Counter()); qlen={}; rlen={}
for l in open("aln.paf"):
    p=l.split("\t"); q,ql,qs,qe,st,r,rl,rs,re_,m,bl,mq=p[:12]
    ql=int(ql); qlen[q]=ql; rlen[r]=int(rl)
    if int(bl)<50000 or int(mq)<30: continue
    hits[q][r]+=int(qe)-int(qs)
rows=[]
for q in sorted(qlen,key=lambda x:-qlen[x]):
    if qlen[q]<1e6: continue
    c=hits[q]; tot=sum(c.values()); best,bp=(c.most_common(1)[0] if c else ("-",0))
    second=c.most_common(2)[1] if len(c)>1 else ("-",0)
    rows.append((q,qlen[q],best,bp,100*bp/qlen[q] if qlen[q] else 0,second[0],second[1],100*second[1]/qlen[q]))
with open("scaffold_to_GopEvg.tsv","w") as o:
    o.write("scaffold\tlength\tbest_ref_chr\tbest_bp\tbest_pct_of_scaffold\tsecond_ref_chr\tsecond_bp\tsecond_pct\n")
    for r in rows: o.write("\t".join(f"{x:.1f}" if isinstance(x,float) else str(x) for x in r)+"\n")
for r in rows[:30]: print(f"{r[0]:14s} {r[1]/1e6:7.1f} Mb  -> {r[2]:14s} {r[4]:5.1f}%   second {r[5]:14s} {r[7]:5.1f}%")
# which reference chromosomes are covered by >1 scaffold (fusion/split check)
inv=collections.defaultdict(list)
for r in rows:
    if r[4]>20: inv[r[2]].append((r[0],r[4]))
print("\nreference chromosomes hit by >1 scaffold (>20% of scaffold):")
for k,v in sorted(inv.items()):
    if len(v)>1: print(" ",k,rlen.get(k,0)/1e6,"Mb:",v)
PY
# dot plot
python3 - <<PY
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, collections
qlen={}; rlen={}; segs=[]
for l in open("aln.paf"):
    p=l.split("\t"); 
    if int(p[10])<20000 or int(p[11])<30: continue
    qlen[p[0]]=int(p[1]); rlen[p[5]]=int(p[6]); segs.append((p[0],int(p[2]),int(p[3]),p[5],int(p[7]),int(p[8]),p[4]))
Q=[q for q in sorted(qlen,key=lambda x:-qlen[x]) if qlen[q]>=5e6]; R=sorted(rlen,key=lambda x:-rlen[x])
qo={}; x=0
for q in Q: qo[q]=x; x+=qlen[q]
ro={}; y=0
for r in R: ro[r]=y; y+=rlen[r]
fig,ax=plt.subplots(figsize=(11,11))
for q,qs,qe,r,rs,re_,st in segs:
    if q not in qo: continue
    ys=(rs,re_) if st=="+" else (re_,rs)
    ax.plot([qo[q]+qs,qo[q]+qe],[ro[r]+ys[0],ro[r]+ys[1]],color="#2E3F6E" if st=="+" else "#B5791F",lw=0.8)
for q in Q: ax.axvline(qo[q],color="#ccc",lw=0.4)
for r in R: ax.axhline(ro[r],color="#ccc",lw=0.4)
ax.set_xticks([qo[q]+qlen[q]/2 for q in Q]); ax.set_xticklabels([q.replace("scaffold_","s") for q in Q],rotation=90,fontsize=6)
ax.set_yticks([ro[r]+rlen[r]/2 for r in R]); ax.set_yticklabels([r.replace("GopEvg_","") for r in R],fontsize=6)
ax.set_xlim(0,x); ax.set_ylim(0,y); ax.set_xlabel("M. impressa v3 ($LABEL) scaffolds >= 5 Mb"); ax.set_ylabel("G. evgoodei rGopEvg1 chromosomes")
ax.set_title("Whole-genome alignment (minimap2 asm20; blue +, amber −)")
plt.tight_layout(); plt.savefig("dotplot_$LABEL.png",dpi=200)
PY
echo "== DONE $(date)"
