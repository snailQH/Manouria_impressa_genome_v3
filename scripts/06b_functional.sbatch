#!/bin/bash
#SBATCH --job-name=functional
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06b_functional_%A.out
# Step 14: annotation QC + functional annotation of the BRAKER3 gene set.
#   gene-set statistics (python); BUSCO v5.8.0 protein mode vertebrata_odb10; eggNOG-mapper (DIAMOND, eggNOG 5.0.2: COG/GO/KEGG/Pfam/preferred name).
#   usage: sbatch 06b_functional.sbatch <label>
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; L=$1; A=$NV/05_annotation/$L; B=$A/02_braker3; [ -s $B/braker.gff3 ] || B=$A/02_braker3/braker; O=$A/03_functional; mkdir -p $O; cd $O
echo "== $L $(date) node=$(hostname)"; ls -la $B/braker.gff3 $B/braker.aa $B/braker.codingseq
# longest isoform per gene for BUSCO/eggNOG
python3 - $B/braker.aa $B/braker.gff3 <<'PY'
import sys,re,collections
aa,gff=sys.argv[1:3]; seqs={}; cur=None
for l in open(aa):
    if l.startswith(">"): cur=l[1:].split()[0]; seqs[cur]=[]
    else: seqs[cur].append(l.strip())
seqs={k:"".join(v).rstrip("*") for k,v in seqs.items()}
genes=collections.defaultdict(list)
for k in seqs: genes[k.rsplit(".",1)[0]].append(k)
with open("proteins.longest.faa","w") as f:
    for g,ts in genes.items():
        t=max(ts,key=lambda x:len(seqs[x])); f.write(f">{t}\n{seqs[t]}\n")
ex=collections.Counter(); glen=[]; cds=collections.defaultdict(int); intr=collections.defaultdict(list); last={}
for l in open(gff):
    if l.startswith("#"): continue
    p=l.rstrip("\n").split("\t");
    if len(p)<9: continue
    if p[2]=="gene": glen.append(int(p[4])-int(p[3])+1)
    if p[2]=="exon":
        t=re.search(r"Parent=([^;]+)",p[8]).group(1); ex[t]+=1
        if t in last: intr[t].append(int(p[3])-last[t]-1)
        last[t]=int(p[4])
    if p[2]=="CDS": t=re.search(r"Parent=([^;]+)",p[8]).group(1); cds[t]+=int(p[4])-int(p[3])+1
import statistics as st
ne=list(ex.values()); il=[i for v in intr.values() for i in v if i>0]; cl=list(cds.values())
with open("gene_set_stats.tsv","w") as f:
    for k,v in [("genes",len(genes)),("transcripts",len(seqs)),("mean_transcripts_per_gene",f"{len(seqs)/len(genes):.2f}"),
        ("mean_gene_length_bp",f"{st.mean(glen):.0f}"),("mean_CDS_length_bp",f"{st.mean(cl):.0f}"),("mean_exons_per_transcript",f"{st.mean(ne):.2f}"),
        ("mean_exon_length_bp",f"{st.mean(cl)/st.mean(ne):.0f}"),("mean_intron_length_bp",f"{st.mean(il):.0f}"),("median_protein_length_aa",st.median(len(s) for s in seqs.values())),
        ("single_exon_transcripts",sum(1 for v in ne if v==1))]: f.write(f"{k}\t{v}\n")
print(open("gene_set_stats.tsv").read())
PY
conda activate /data1/opt/conda/envs/busco
busco -i proteins.longest.faa -l vertebrata_odb10 -m proteins -c 48 -o busco_prot_vertebrata --download_path /data1/reference/busco_downloads/busco_downloads --offline > busco_prot.log 2>&1 || true
grep -E "C:|Complete|Fragmented|Missing" busco_prot_vertebrata/short_summary*.txt | head -5
conda activate /data1/opt/conda/envs/eggnog
emapper.py -i proteins.longest.faa -o eggnog --cpu 48 --data_dir /data1/reference/eggnog -m diamond --dmnd_db /data1/reference/eggnog/eggnog_proteins.dmnd --tax_scope 7742 --go_evidence non-electronic --override > emapper.log 2>&1 || true
N=$(grep -vc "^#" eggnog.emapper.annotations 2>/dev/null || echo 0); T=$(grep -c ">" proteins.longest.faa)
echo "eggNOG: $N of $T genes annotated"; awk -F'\t' '!/^#/{if($10!="-")go++; if($12!="-")ko++; if($21!="-")pf++; n++} END{print "GO:",go," KEGG_ko:",ko," Pfam:",pf," total:",n}' eggnog.emapper.annotations 2>/dev/null
echo "== DONE $(date)"
