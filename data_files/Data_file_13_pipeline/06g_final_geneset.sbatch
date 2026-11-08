#!/bin/bash
#SBATCH --job-name=final_geneset
#SBATCH --partition=main
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/06g_final_geneset_%A.out
# Step 14f: final gene set = AUGUSTUS-with-hints genes (26,279; BUSCO-protein C 97.6) that carry at least one line of evidence:
#   RNA-seq read coverage (>= 1 heart/lung read over the gene span) OR an eggNOG-mapper ortholog hit OR a complete/fragmented BUSCO.
#   Pre-registered before the numbers were seen. Outputs: Mimp_v3.genes.gtf/gff3, .proteins.faa, .cds.fna (release names, IDs MIMP3_gNNNNN),
#   gene-set statistics, BUSCO-protein, eggNOG annotation rate, RNA-seq support, Table 3 inputs.
set -euo pipefail
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; A=$NV/05_annotation/Mimp_v3; B=$A/02_braker3; C=$A/05_annotation_check; O=$A/06_final_geneset; mkdir -p $O; cd $O
G=$A/01_repeats/Mimp_v3.softmasked.fa; SIF=/data1/reference/containers/braker3.sif; RUN="singularity exec -B /data2,/data1 $SIF"
source /data1/opt/conda/etc/profile.d/conda.sh
# 1. eggNOG-mapper on the AUGUSTUS set (longest isoform per gene)
conda activate /data1/opt/conda/envs/eggnog
[ -s eggnog_aug.emapper.annotations ] || emapper.py -i $C/augustus.longest.faa -o eggnog_aug --cpu 48 --data_dir /data1/reference/eggnog -m diamond --dmnd_db /data1/reference/eggnog/eggnog_proteins.dmnd --tax_scope 7742 --go_evidence non-electronic --override > emapper_aug.log 2>&1
grep -v "^#" eggnog_aug.emapper.annotations | cut -f1 | sed 's/\.t[0-9]*$//' | sort -u > genes_with_eggnog.txt
# 2. RNA-seq support (>= 1 read) from the decision job; BUSCO hits
conda activate /data1/opt/conda/envs/assembly
awk '$7+$8>=1{print $4}' $C/chosen_genes_rnaseq_counts.bed | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort -u > genes_with_rnaseq.txt
cat $C/busco_augustus/run_vertebrata_odb10/full_table.tsv | awk -F'\t' '$2=="Complete"||$2=="Duplicated"||$2=="Fragmented"{print $3}' | sed 's/\.t[0-9]*$//' | sort -u > genes_with_busco.txt
cat genes_with_eggnog.txt genes_with_rnaseq.txt genes_with_busco.txt | sort -u > genes_supported.txt
awk -F'\t' '$3=="gene"{print $9}' $B/Augustus/augustus.hints.gtf | sort -u > genes_all.txt   # AUGUSTUS gene lines carry the bare id
echo "AUGUSTUS genes: $(wc -l < genes_all.txt); with eggNOG hit: $(wc -l < genes_with_eggnog.txt); with RNA-seq: $(wc -l < genes_with_rnaseq.txt); with BUSCO: $(wc -l < genes_with_busco.txt); supported (any): $(wc -l < genes_supported.txt); unsupported (dropped): $(comm -23 genes_all.txt genes_supported.txt | wc -l)"
# 3. filtered GTF → renamed GFF3 / proteins / CDS
python3 - <<'PY'
import re, collections
keep=set(l.strip() for l in open('genes_supported.txt'))
gtf='/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/05_annotation/Mimp_v3/02_braker3/Augustus/augustus.hints.gtf'
# order genes by chromosome then position for stable IDs
genes=collections.OrderedDict(); lines=[]
for l in open(gtf):
    if l.startswith('#'): continue
    p=l.rstrip('\n').split('\t')
    if len(p)<9: continue
    m=re.search(r'gene_id "([^"]+)"',p[8])
    g=m.group(1) if m else (p[8].strip() if p[2]=='gene' else p[8].strip().rsplit('.t',1)[0])   # gene: 'g1'; transcript: 'g1.t1'; others: attributes
    if g not in keep: continue
    if p[2]=='gene': genes[g]=(p[0],int(p[3]))
    lines.append(p)
def key(c):
    m=re.match(r'chr(\d+)$',c); return (0,int(m.group(1))) if m else (1,int(re.search(r'(\d+)$',c).group(1)))
order=sorted(genes,key=lambda g:(key(genes[g][0]),genes[g][1])); newid={g:f'MIMP3_g{i:05d}' for i,g in enumerate(order,1)}
with open('Mimp_v3.genes.gtf','w') as o:
    for p in lines:
        a=p[8]
        for g,n in newid.items(): pass
        if p[2]=='gene': a=newid[a.strip()]
        elif p[2]=='transcript': g,t=a.strip().rsplit('.t',1); a=f'{newid[g]}.t{t}'
        else:
            a=re.sub(r'gene_id "([^"]+)"',lambda m:f'gene_id "{newid[m.group(1)]}"',a)
            a=re.sub(r'transcript_id "([^"]+)\.t(\d+)"',lambda m:f'transcript_id "{newid[m.group(1)]}.t{m.group(2)}"',a)
        o.write('\t'.join(p[:8]+[a])+'\n')
with open('gene_id_map.tsv','w') as o:
    o.write('release_id\taugustus_id\tchromosome\tstart\n'); [o.write(f'{newid[g]}\t{g}\t{genes[g][0]}\t{genes[g][1]}\n') for g in order]
print('written genes:',len(order))
PY
# proteins / CDS from the GTF (gffread uses the CDS lines); GFF3 written directly so that gene and mRNA features are present (GWH needs them)
$RUN /opt/ETP/tools/gffread Mimp_v3.genes.gtf -g $G -y Mimp_v3.proteins.raw.faa -x Mimp_v3.cds.fna 2> gffread_final.log || true
awk '/^>/{print; next}{gsub(/[.*]/,""); print}' Mimp_v3.proteins.raw.faa > Mimp_v3.proteins.faa   # gffread writes '.' for frame gaps and '*' at stops; neither is accepted by DIAMOND/GWH
python3 - <<'PY2'
import re, collections
feat=collections.defaultdict(list); genes={}; trans={}
for l in open('Mimp_v3.genes.gtf'):
    p=l.rstrip('\n').split('\t')
    if len(p)<9: continue
    if p[2]=='gene': genes[p[8].strip()]=(p[0],int(p[3]),int(p[4]),p[6])
    elif p[2]=='transcript': t=p[8].strip(); trans[t]=(p[0],int(p[3]),int(p[4]),p[6],t.rsplit('.t',1)[0])
    elif p[2] in ('exon','CDS','start_codon','stop_codon','intron'):
        t=re.search(r'transcript_id "([^"]+)"',p[8]).group(1); feat[t].append((p[2],int(p[3]),int(p[4]),p[6],p[7]))
with open('Mimp_v3.genes.gff3','w') as o:
    o.write('##gff-version 3\n##genome Mimp_v3 (Manouria impressa); source BRAKER 3.0.8 / AUGUSTUS with RNA-seq + OrthoDB hints; evidence-supported gene set\n')
    for g,(c,s_,e,st) in genes.items():
        o.write(f'{c}\tAUGUSTUS\tgene\t{s_}\t{e}\t.\t{st}\t.\tID={g};Name={g}\n')
        for t,(tc,ts,te,tst,tg) in trans.items():
            if tg!=g: continue
            o.write(f'{tc}\tAUGUSTUS\tmRNA\t{ts}\t{te}\t.\t{tst}\t.\tID={t};Parent={g}\n')
            n=0
            for typ,fs,fe,fst,ph in sorted(feat[t],key=lambda x:x[1]):
                if typ in ('exon','CDS'):
                    n+=1 if typ=='exon' else 0; fid=f'{t}.exon{n}' if typ=='exon' else f'{t}.cds'
                    o.write(f'{tc}\tAUGUSTUS\t{typ}\t{fs}\t{fe}\t.\t{fst}\t{ph if typ=="CDS" else "."}\tID={fid};Parent={t}\n')
print('gff3: genes',len(genes),'mRNAs',len(trans))
PY2
# gff3 header + sanity: features within sequence bounds, IDs unique
python3 - <<'PY'
L={l.split()[0]:int(l.split()[1]) for l in open('/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/06_release/Mimp_v3.fa.fai')}
bad=0; ids=set(); dup=0; n=0
for l in open('Mimp_v3.genes.gff3'):
    if l.startswith('#'): continue
    p=l.rstrip('\n').split('\t'); n+=1
    if p[0] not in L or int(p[4])>L[p[0]] or int(p[3])<1: bad+=1
    if p[2] in ('gene','mRNA'):
        i=p[8].split('ID=')[1].split(';')[0]
        if i in ids: dup+=1
        ids.add(i)
print(f'gff3 features {n}; out of bounds {bad}; duplicate IDs {dup}')
PY
# 4. statistics, BUSCO, eggNOG rate, RNA-seq support on the final set
python3 - <<'PY'
import re, collections, statistics as st
seqs={}; cur=None
for l in open('Mimp_v3.proteins.faa'):
    if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
    else: seqs[cur].append(l.strip())
seqs={k:''.join(v).rstrip('.').rstrip('*') for k,v in seqs.items()}
genes=collections.defaultdict(list)
for k in seqs: genes[k.rsplit('.t',1)[0]].append(k)
with open('proteins.longest.faa','w') as f:
    for g,ts in genes.items(): t=max(ts,key=lambda x:len(seqs[x])); f.write(f'>{t}\n{seqs[t]}\n')
ex=collections.Counter(); cds=collections.defaultdict(int); glen=[]; intr=[]; last={}
for l in open('Mimp_v3.genes.gff3'):
    if l.startswith('#'): continue
    p=l.rstrip('\n').split('\t')
    if p[2]=='gene': glen.append(int(p[4])-int(p[3])+1)
    if p[2]=='exon':
        t=re.search(r'Parent=([^;]+)',p[8]).group(1); ex[t]+=1
        if t in last: intr.append(int(p[3])-last[t]-1)
        last[t]=int(p[4])
    if p[2]=='CDS': t=re.search(r'Parent=([^;]+)',p[8]).group(1); cds[t]+=int(p[4])-int(p[3])+1
ne=list(ex.values()); il=[i for i in intr if i>0]; cl=list(cds.values())
rows=[('genes',len(genes)),('transcripts',len(seqs)),('mean_transcripts_per_gene',f'{len(seqs)/len(genes):.2f}'),('mean_gene_length_bp',f'{st.mean(glen):.0f}'),('mean_CDS_length_bp',f'{st.mean(cl):.0f}'),('mean_exons_per_transcript',f'{st.mean(ne):.2f}'),('mean_exon_length_bp',f'{st.mean(cl)/st.mean(ne):.0f}'),('mean_intron_length_bp',f'{st.mean(il):.0f}'),('median_protein_length_aa',st.median(len(s) for s in seqs.values())),('single_exon_transcripts',sum(1 for v in ne if v==1))]
open('gene_set_stats.tsv','w').write('\n'.join(f'{k}\t{v}' for k,v in rows)+'\n'); print(open('gene_set_stats.tsv').read())
PY
conda activate /data1/opt/conda/envs/busco
busco -i proteins.longest.faa -l vertebrata_odb10 -m proteins -c 48 -o busco_final --download_path /data1/reference/busco_downloads/busco_downloads --offline > busco_final.log 2>&1 || true
echo "== final BUSCO-protein: $(grep -h 'C:' busco_final/short_summary*.txt)"
# eggNOG rate on the final set (map old ids → new)
python3 - <<'PY'
m={l.split('\t')[1]:l.split('\t')[0] for l in open('gene_id_map.tsv') if not l.startswith('release_id')}
n=0; got=set(); go=ko=pf=0
for l in open('eggnog_aug.emapper.annotations'):
    if l.startswith('#'): continue
    p=l.rstrip('\n').split('\t'); g=p[0].rsplit('.t',1)[0]
    if g in m and g not in got:
        got.add(g); go+=p[9]!='-'; ko+=p[11]!='-'; pf+=p[20]!='-'
print(f'eggNOG on final set: {len(got)} of {len(m)} genes ({100*len(got)/len(m):.1f} %); GO {go}; KEGG ko {ko}; Pfam {pf}')
open('eggnog_rate.txt','w').write(f'{len(got)}\t{len(m)}\t{go}\t{ko}\t{pf}\n')
PY
conda activate /data1/opt/conda/envs/assembly
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{match($9,/ID=[^;]+/); print $1,$4-1,$5,substr($9,RSTART+3,RLENGTH-3),".",$7}' Mimp_v3.genes.gff3 | sort -k1,1 -k2,2n > final_genes.bed
bedtools multicov -bams $A/rnaseq/IT3.bam $A/rnaseq/IT4.bam -bed final_genes.bed > final_genes_rnaseq_counts.bed
awk '{t++; if($7+$8>=10) s++; if($7+$8>=1) a++} END{printf "RNA-seq support (final): >=10 reads %d/%d (%.1f %%); >=1 read %d (%.1f %%)\n", s,t,100*s/t,a,100*a/t}' final_genes_rnaseq_counts.bed | tee rnaseq_support_final.txt
awk -F'\t' '$3=="gene"{c[$1 ~ /^chr/]++} END{print "genes on chromosomes:",c[1],"unplaced:",c[0]+0}' Mimp_v3.genes.gff3
cp Mimp_v3.genes.gff3 Mimp_v3.proteins.faa Mimp_v3.cds.fna gene_id_map.tsv $NV/06_release/; md5sum $NV/06_release/Mimp_v3.genes.gff3 $NV/06_release/Mimp_v3.proteins.faa $NV/06_release/Mimp_v3.cds.fna >> $NV/06_release/release_md5.txt
cp gene_set_stats.tsv eggnog_rate.txt rnaseq_support_final.txt $NV/09_manuscript_inputs/; cp busco_final/short_summary*.txt $NV/09_manuscript_inputs/busco_protein_final.txt; cp eggnog_aug.emapper.annotations $NV/06_release/Mimp_v3.eggnog.annotations.tsv
echo "== DONE $(date)"
