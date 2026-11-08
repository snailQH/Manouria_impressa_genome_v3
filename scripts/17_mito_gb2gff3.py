#!/usr/bin/env python3
"""Convert the MitoHiFi GenBank annotation of the mitogenome to GFF3 for the GWH organelle deposit (and Data file 7).
Pure Python (no Biopython on the cluster envs). Features: gene -> CDS (transl_table=2) | tRNA | rRNA; strand from complement();
join() locations become multi-line features sharing one ID. Repeated names (the tRNA-Thr2/-Pro2 copies of the tandem array) get a _2/_3 suffix. The seqid is taken from the release FASTA header, not from the GenBank LOCUS.
usage: 17_mito_gb2gff3.py Mimp_v3_mito.gb Mimp_v3_mito.fa > Mimp_v3_mito.gff3"""
import re, sys
gb, fa = sys.argv[1], sys.argv[2]
seqid = open(fa).readline()[1:].split()[0]
seqlen = sum(len(l.strip()) for l in open(fa) if not l.startswith('>'))
lines = open(gb).read().split('\n')
feats = []; cur = None
for l in lines:
    if l.startswith('ORIGIN'): break
    m = re.match(r'^ {5}(\S+) {2,}(\S.*)$', l)
    if m:
        cur = {'type': m.group(1), 'loc': m.group(2).strip(), 'q': {}}; feats.append(cur); continue
    if cur is not None and l.startswith(' ' * 21):
        t = l.strip()
        if t.startswith('/'):
            k, _, v = t[1:].partition('='); cur['q'][k] = v.strip('"'); cur['last'] = k
        elif 'last' in cur and cur['q'][cur['last']] is not None:
            cur['q'][cur['last']] += t.rstrip('"')
        elif not t.startswith('/') and 'last' not in cur:   # continuation of a long location
            cur['loc'] += t
def parse_loc(loc):
    strand = '+'
    if loc.startswith('complement('): strand = '-'; loc = loc[len('complement('):-1]
    if loc.startswith('join('): loc = loc[5:-1]
    parts = []
    for p in loc.split(','):
        p = p.strip().replace('<', '').replace('>', '')
        a, b = (p.split('..') + [p])[:2]; parts.append((int(a), int(b)))
    return strand, parts
out = ['##gff-version 3', f'##sequence-region {seqid} 1 {seqlen}', f'{seqid}\tMitoHiFi\tregion\t1\t{seqlen}\t.\t+\t.\tID={seqid};Name={seqid};Is_circular=true;genome=mitochondrion;mol_type=genomic DNA;organism=Manouria impressa']
genes = [f for f in feats if f['type'] == 'gene']; children = [f for f in feats if f['type'] in ('CDS', 'tRNA', 'rRNA')]
names = {}
def name_of(f):
    n = f['q'].get('gene') or f['q'].get('product') or f['type']
    return n
for i, g in enumerate(genes, 1):
    strand, parts = parse_loc(g['loc']); gid = f'MIMP3_mt_g{i:02d}'; n = name_of(g)
    names.setdefault(n, 0); names[n] += 1
    child = next((c for c in children if c['loc'] == g['loc'] and c not in [x for x in genes]), None)
    if child is None: child = next((c for c in children if parse_loc(c['loc'])[1][0][0] == parts[0][0]), None)
    label = n if names[n] == 1 else f'{n}_{names[n]}'          # MitoHiFi already names second isoacceptors Leu2/Ser2; array repeats get _2, _3 …
    note = ';Note=unit of the tandem tRNA-Thr/tRNA-Pro array adjacent to the control region (four units in the deposited consensus)' if n in ('tRNA-Thr', 'tRNA-Pro', 'tRNA-Thr2', 'tRNA-Pro2') and parts[0][0] > 15000 else ''
    s, e = min(p[0] for p in parts), max(p[1] for p in parts)
    biotype = {'CDS': 'protein_coding', 'tRNA': 'tRNA', 'rRNA': 'rRNA'}[child['type']]
    out.append(f'{seqid}\tMitoHiFi\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid};Name={label};gene={n};gene_biotype={biotype}{note}')
    cid = gid.replace('_g', '_' + {'CDS': 'cds', 'tRNA': 't', 'rRNA': 'r'}[child['type']])
    prod = child['q'].get('product') or child['q'].get('gene') or n
    extra = ';transl_table=2' if child['type'] == 'CDS' else ''
    phase_pos = 0
    for (a, b) in (parts if strand == '+' else parts[::-1]):
        phase = str((3 - phase_pos % 3) % 3) if child['type'] == 'CDS' else '.'
        out.append(f'{seqid}\tMitoHiFi\t{child["type"]}\t{a}\t{b}\t.\t{strand}\t{phase}\tID={cid};Parent={gid};Name={label};product={prod};gene={n}{extra}')
        phase_pos += b - a + 1
    if child['type'] == 'CDS' and child['q'].get('translation', '').rstrip('*') and len(parts) == 1:
        aa = child['q']['translation'].rstrip('*'); nt = e - s + 1
        if abs(nt - 3 * len(aa)) > 3: print(f'# warning {n}: {nt} nt vs {len(aa)} aa', file=sys.stderr)
# control region: the unannotated span between the last tRNA (Pro) and tRNA-Phe at position 1, as in all vertebrate mitogenomes
last_end = max(max(b for a, b in parse_loc(g['loc'])[1]) for g in genes)
if seqlen - last_end > 500:
    out.append(f'{seqid}\tMitoHiFi\tD_loop\t{last_end + 1}\t{seqlen}\t.\t+\t.\tID=MIMP3_mt_dloop;Name=control region;Note=inferred from position between tRNA-Pro and tRNA-Phe, not annotated by MitoHiFi')
print('\n'.join(out))
