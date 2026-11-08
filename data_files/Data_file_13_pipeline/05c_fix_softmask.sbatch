#!/bin/bash
#SBATCH --job-name=softmask_fix
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/05c_softmask_fix_%A.out
# Step 12c: correct soft-masked genome = union of both RepeatMasker passes (the second pass uppercases its input, so
#   `cp rm_dfam/*.masked.masked` kept only the Dfam masking), restricted to the FCS-cleaned release and renamed to release names.
#   Output: 05_annotation/Mimp_v3/01_repeats/Mimp_v3.softmasked.fa  (the genome BRAKER3 annotates).
set -euo pipefail
source /data1/opt/conda/etc/profile.d/conda.sh; conda activate /data1/opt/conda/envs/assembly
NV=/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version; R=$NV/05_annotation/Mimp_v3/01_repeats; cd $R
mv genome.softmasked.fa genome.softmasked.dfam-only.WRONG.fa 2>/dev/null || true
# union: species-library soft-masked FASTA + Dfam-pass intervals (GFF, pre-FCS names)
grep -v "^#" rm_dfam/genome.fa.masked.out.gff | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' > rm_dfam/dfam_pass.bed
bedtools maskfasta -soft -fi rm_custom/genome.fa.masked -bed rm_dfam/dfam_pass.bed -fo genome.softmasked.pre-fcs.fa
# restrict to the cleaned release, rename unplaced scaffolds to release names
python3 - <<'PY'
NV='/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version'
m={}
for l in open(f'{NV}/06_release/Mimp_v3.names.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0]!='release_name': m[p[-1]]=p[0]        # pre_fcs_name -> release_name
out=open('Mimp_v3.softmasked.fa','w'); keep=False; n=0; low=0; tot=0
for l in open('genome.softmasked.pre-fcs.fa'):
    if l.startswith('>'):
        name=l[1:].split()[0]; keep=name in m
        if keep: out.write('>'+m[name]+'\n'); n+=1
    elif keep:
        out.write(l); s=l.strip(); tot+=len(s); low+=sum(1 for c in s if c.islower())
out.close(); print(f'Mimp_v3.softmasked.fa: {n} sequences, {tot:,} bp, soft-masked {low:,} bp = {100*low/tot:.2f} %')
PY
samtools faidx Mimp_v3.softmasked.fa; cmp <(cut -f1,2 Mimp_v3.softmasked.fa.fai | sort) <(cut -f1,2 $NV/06_release/Mimp_v3.fa.fai | sort) && echo "names/lengths identical to 06_release/Mimp_v3.fa"
# repeats GFF for the release (renamed, contaminants dropped)
python3 - <<'PY'
NV='/data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version'; m={}
for l in open(f'{NV}/06_release/Mimp_v3.names.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0]!='release_name': m[p[-1]]=p[0]
n=0
with open('Mimp_v3.repeats.gff','w') as o:
    for l in open('repeats.gff'):
        p=l.split('\t')
        if p[0] in m: p[0]=m[p[0]]; o.write('\t'.join(p)); n+=1
print('repeat features on the release:', n)
PY
echo "== DONE $(date)"
