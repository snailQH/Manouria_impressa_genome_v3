#!/bin/bash
#SBATCH --job-name=fcs_clean
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/13_fcs_clean_%A.out
# Step 11f: apply the FCS-GX action reports (EXCLUDE / TRIM / FIX) to the three release assemblies with `gx clean-genome`,
#   renumber unplaced / haplotype scaffolds by length, rebuild AGP + names map + haplotype→primary tables + stats + md5.
#   The pre-screen files are kept as *.pre-fcs.* for provenance. Contaminant sequences are written to 07_contamination/<label>/gx/contam.fa.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; C=/data1/reference/containers; R=$NV/06_release; cd $R
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
for L in Mimp_v3 Mimp_v3_hap1 Mimp_v3_hap2; do
  REP=$(ls $NV/07_contamination/$L/gx/*.fcs_gx_report.txt)
  [ -s $L.pre-fcs.fa ] || cp $L.fa $L.pre-fcs.fa
  singularity exec -B /data2,/data1 $C/fcs-gx.sif gx clean-genome -i $L.pre-fcs.fa --action-report $REP --contam-fasta-out $NV/07_contamination/$L/gx/contam.fa -o $L.cleaned.tmp.fa
  echo "$L: $(grep -c '>' $L.pre-fcs.fa) → $(grep -c '>' $L.cleaned.tmp.fa) sequences; contaminant: $(grep -c '>' $NV/07_contamination/$L/gx/contam.fa) sequences, $(grep -v '>' $NV/07_contamination/$L/gx/contam.fa | tr -d '\n' | wc -c) bp"
done
python3 - <<'PY'
import subprocess, re, collections
R='/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/06_release'
def read(f):
    seqs=[]; name=None; buf=[]
    for l in open(f):
        if l.startswith('>'):
            if name: seqs.append((name,''.join(buf)))
            name=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if name: seqs.append((name,''.join(buf)))
    return seqs
def write(f, seqs):
    with open(f,'w') as o:
        for n,s in seqs:
            o.write('>'+n+'\n'); [o.write(s[i:i+100]+'\n') for i in range(0,len(s),100)]
# ---- primary: chromosomes keep names; unplaced renumbered by length; split pieces (from FIX) get a suffix in the map
seqs=read(f'{R}/Mimp_v3.cleaned.tmp.fa'); chrs=[(n,s) for n,s in seqs if n.startswith('chr')]
unp=sorted([(n,s) for n,s in seqs if not n.startswith('chr')], key=lambda x:-len(x[1]))
newmap={}; out=list(chrs)
for i,(n,s) in enumerate(unp,1): newmap[n]=f'scaffold_{i}'; out.append((f'scaffold_{i}',s))
write(f'{R}/Mimp_v3.fa', out)
old=dict((l.split('\t')[0],l.rstrip('\n').split('\t')) for l in open(f'{R}/Mimp_v3.names.tsv'))   # release_name -> row (curated_name,...)
hdr=open(f'{R}/Mimp_v3.names.tsv').readline().rstrip('\n')
with open(f'{R}/Mimp_v3.names.tsv','w') as o:
    o.write(hdr+'\tpre_fcs_name\n')
    for n,s in out:
        src=n if n.startswith('chr') else [k for k,v in newmap.items() if v==n][0]
        row=old.get(src.split('_')[0]+'_'+src.split('_')[1] if '_' in src else src, [src,'-','-','-','-','-'])
        o.write('\t'.join([n]+row[1:])+f'\t{src}\n')
# AGP: keep chromosome objects; rewrite unplaced objects with new names; drop removed
agp=[l for l in open(f'{R}/Mimp_v3.pre-fcs.agp')] if __import__('os').path.exists(f'{R}/Mimp_v3.pre-fcs.agp') else [l for l in open(f'{R}/Mimp_v3.agp')]
open(f'{R}/Mimp_v3.pre-fcs.agp','w').writelines(agp)
keep=set(n for n,_ in out)
with open(f'{R}/Mimp_v3.agp','w') as o:
    for l in agp:
        if l.startswith('#'): o.write(l); continue
        p=l.rstrip('\n').split('\t'); obj=p[0]
        if obj.startswith('chr'): o.write(l)
        elif obj in newmap and len([1 for n,s in unp if n==obj])==1: p[0]=newmap[obj]; o.write('\t'.join(p)+'\n')
print('primary written:', len(out), 'sequences;', len(chrs), 'chromosomes;', len(unp), 'unplaced')
# ---- haplotypes: renumber by length
for h in ('1','2'):
    seqs=read(f'{R}/Mimp_v3_hap{h}.cleaned.tmp.fa'); seqs=sorted(seqs, key=lambda x:-len(x[1]))
    with open(f'{R}/Mimp_v3_hap{h}.names.tsv','w') as o:
        o.write('release_name\tpre_fcs_name\tlength\n')
        out=[]
        for i,(n,s) in enumerate(seqs,1): out.append((f'h{h}_scaffold_{i}',s)); o.write(f'h{h}_scaffold_{i}\t{n}\t{len(s)}\n')
    write(f'{R}/Mimp_v3_hap{h}.fa', out); print(f'hap{h} written:', len(out))
PY
rm -f *.cleaned.tmp.fa
for L in Mimp_v3 Mimp_v3_hap1 Mimp_v3_hap2; do samtools faidx $L.fa; python3 $NV/scripts/asm_stats.py $L.fa > $L.stats.tsv; echo "== $L"; head -n 8 $L.stats.tsv; done
# haplotype → primary tables against the cleaned primary
for H in 1 2; do OUT=Mimp_v3_hap$H; minimap2 -x asm5 -t 16 --secondary=no Mimp_v3.fa $OUT.fa > ${OUT}_to_primary.paf 2>/dev/null
python3 - $OUT <<'PY'
import sys, collections
o=sys.argv[1]; hits=collections.defaultdict(collections.Counter); L={}
for l in open(o+"_to_primary.paf"):
    p=l.split("\t"); L[p[0]]=int(p[1])
    if int(p[10])<20000 or int(p[11])<30: continue
    hits[p[0]][p[5]]+=int(p[3])-int(p[2])
with open(o+"_to_primary.tsv","w") as f:
    f.write("hap_scaffold\tlength\tbest_primary_sequence\taligned_bp\tpct_of_scaffold\n")
    for s in sorted(L,key=lambda x:-L[x]):
        b=hits[s].most_common(1)[0] if hits[s] else ("-",0); f.write(f"{s}\t{L[s]}\t{b[0]}\t{b[1]}\t{100*b[1]/L[s]:.1f}\n")
PY
done
# telomeres unchanged (chromosomes only); re-derive the contaminant list table for the supplement
(echo -e "assembly\tsequence_pre_fcs\tlength_bp\taction\tdivision\tcoverage\ttop_taxon"; for L in Mimp_v3 Mimp_v3_hap1 Mimp_v3_hap2; do grep -v "^#" $NV/07_contamination/$L/gx/*.fcs_gx_report.txt | awk -F'\t' -v L=$L 'BEGIN{OFS="\t"}{print L,$1,$3-$2+1,$5,$6,$7,$8}'; done) > fcs_gx_actions_applied.tsv
md5sum Mimp_v3.fa Mimp_v3.agp Mimp_v3_hap1.fa Mimp_v3_hap2.fa > release_md5.txt; cat release_md5.txt
# RepeatMasker summary restricted to the kept primary sequences (contaminant scaffolds dropped from the .out)
RM=$NV/05_annotation/Mimp_v3/01_repeats; cut -f1 Mimp_v3.names.tsv | tail -n +2 > /dev/null
python3 - <<'PY'
R='/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/06_release'; keep={}
for l in open(f'{R}/Mimp_v3.names.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0]!='release_name': keep[p[-1]]=p[0]
open(f'{R}/kept_pre_fcs_names.txt','w').write('\n'.join(keep)+'\n'); print('kept', len(keep))
PY
for d in rm_custom rm_dfam; do f=$(ls $RM/$d/*.out | head -1); awk 'NR==FNR{k[$1]=1; next} FNR<=3 || ($5 in k)' kept_pre_fcs_names.txt $f > $RM/$d/genome.cleaned.out; done
G=$(awk '{s+=$2} END{print s}' Mimp_v3.fa.fai)
singularity exec -B /data2,/data1 $C/tetools.sif buildSummary.pl -genome_size $G $RM/rm_custom/genome.cleaned.out > $RM/rm_custom/genome.cleaned.tbl 2>/dev/null || true
grep -E "bases masked|Total interspersed" $RM/rm_custom/genome.cleaned.tbl || true
echo "== DONE $(date)"
