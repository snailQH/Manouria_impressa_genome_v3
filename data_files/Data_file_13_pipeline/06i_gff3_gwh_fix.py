#!/usr/bin/env python3
"""Make the release GFF3 acceptable to the GWH (table2asn-style) validator without changing any protein sequence:
  1. unique IDs for multi-segment CDS features (tx.cds1, tx.cds2 …) — GWH rejects the GFF3-legal shared ID;
  2. partial gene models (AUGUSTUS models cut by a scaffold end or a gap) are flagged NCBI-style:
     partial=true + start_range=.,<start> (left end open) / end_range=<end>,. (right end open) on mRNA and CDS;
  3. the dangling bases that make a partial CDS not a multiple of 3 are trimmed: phase>0 bases at the 5' end (phase -> 0),
     then the remainder at the 3' end — these bases are not translated by AUGUSTUS, so proteins are unchanged;
  4. CDS that begin inside an N-gap are trimmed to the first non-N codon (5' partial); the leading 'X' are removed from the
     protein FASTA and the matching nucleotides from the CDS FASTA.
usage: 06i_gff3_gwh_fix.py Mimp_v3.fa Mimp_v3.genes.gff3 Mimp_v3.proteins.faa Mimp_v3.cds.fna out_prefix
loads the genome into memory (≈ 3 GB RAM). Writes out_prefix.gff3 / .proteins.faa / .cds.fna / .fixlog.tsv"""
import sys, subprocess, collections
fa, gff, faa, fna, out = sys.argv[1:6]
STOP = {'TAA', 'TAG', 'TGA'}
# genome in memory (2.4 GB) — far faster than 293k faidx calls over NFS
GENOME = {}; name = None; buf = []
for l in open(fa):
    if l.startswith('>'):
        if name: GENOME[name] = ''.join(buf).upper()
        name = l[1:].split()[0]; buf = []
    else: buf.append(l.strip())
if name: GENOME[name] = ''.join(buf).upper()
def seq(chrom, s, e): return GENOME[chrom][s - 1:e]
lengths = {k: len(v) for k, v in GENOME.items()}
RC = str.maketrans('ACGTN', 'TGCAN')
lines = [l.rstrip('\n') for l in open(gff)]
feats = []; header = []
for l in lines:
    if l.startswith('#') or not l.strip(): header.append(l); continue
    f = l.split('\t'); a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x); feats.append([f[:8], a])
by_tx = collections.defaultdict(list)
for i, (f, a) in enumerate(feats):
    if f[2] == 'CDS': by_tx[a['Parent']].append(i)
mrna_idx = {a['ID']: i for i, (f, a) in enumerate(feats) if f[2] == 'mRNA'}
log = [['transcript', 'strand', 'issue', 'action']]; trim5 = {}; trim3 = {}
for tx, idxs in by_tx.items():
    idxs.sort(key=lambda i: int(feats[i][0][3])); strand = feats[idxs[0]][0][6]; chrom = feats[idxs[0]][0][0]
    order = idxs if strand == '+' else idxs[::-1]              # transcript 5'->3'
    first, last = order[0], order[-1]
    cds = ''.join(seq(chrom, int(feats[i][0][3]), int(feats[i][0][4])) for i in idxs)
    if strand == '-': cds = cds.translate(RC)[::-1]
    phase = int(feats[first][0][7]) if feats[first][0][7] != '.' else 0
    issues = []; left_open = right_open = False
    # --- 5' end
    cut5 = phase
    if phase > 0: issues.append(f'phase{phase}')
    body = cds[cut5:]
    nN = 0
    while body[nN:nN + 3] and 'N' in body[nN:nN + 3]: nN += 3
    if nN: issues.append(f'Nstart{nN}'); cut5 += nN; body = body[nN:]
    five_partial = bool(issues) or body[:3] != 'ATG'
    if body[:3] != 'ATG' and 'ATG' not in issues: issues.append('nostart')
    at_seq_edge_5 = (strand == '+' and int(feats[first][0][3]) == 1) or (strand == '-' and int(feats[first][0][4]) == lengths[chrom])
    # --- 3' end
    rem = len(body) % 3
    if rem: issues.append(f'len%3={rem}')
    body3 = body[:len(body) - rem] if rem else body
    three_partial = rem > 0 or body3[-3:] not in STOP
    if body3[-3:] not in STOP: issues.append('nostop')
    if not issues: continue
    # apply trims to the terminal CDS segments (coordinates), set phase 0 on the 5' segment
    def shrink(i, n, side):   # side: 'start' of transcript or 'end' of transcript
        f = feats[i][0]
        if (side == 'start') == (strand == '+'): f[3] = str(int(f[3]) + n)
        else: f[4] = str(int(f[4]) - n)
    if cut5: shrink(first, cut5, 'start'); feats[first][0][7] = '0'; trim5[tx] = (cut5, nN)
    if rem: shrink(last, rem, 'end'); trim3[tx] = rem
    # recompute downstream phases from 0
    run = 0
    for i in order:
        feats[i][0][7] = str((3 - run % 3) % 3); run += int(feats[i][0][4]) - int(feats[i][0][3]) + 1
    gmin = min(int(feats[i][0][3]) for i in idxs); gmax = max(int(feats[i][0][4]) for i in idxs)
    left_open = five_partial if strand == '+' else three_partial
    right_open = three_partial if strand == '+' else five_partial
    tags = {'partial': 'true'}
    if left_open: tags['start_range'] = f'.,{gmin}'
    if right_open: tags['end_range'] = f'{gmax},.'
    for i in idxs + [mrna_idx[tx]]: feats[i][1].update(tags)
    if five_partial and at_seq_edge_5: issues.append('at_sequence_edge')
    log.append([tx, strand, ';'.join(issues), f"trim5={cut5}(N={nN}) trim3={rem} partial={'5' if five_partial else ''}{'3' if three_partial else ''}"])
# unique CDS ids
cnt = collections.Counter()
for f, a in feats:
    if f[2] == 'CDS': cnt[a['Parent']] += 1; a['ID'] = f"{a['Parent']}.cds{cnt[a['Parent']]}"
with open(out + '.gff3', 'w') as w:
    w.write('\n'.join(header) + '\n')
    for f, a in feats: w.write('\t'.join(f) + '\t' + ';'.join(f'{k}={v}' for k, v in a.items()) + '\n')
# proteins / CDS fasta: strip leading X / N for the N-start cases
def fasta(p):
    d = collections.OrderedDict(); k = None
    for l in open(p):
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = [l.rstrip('\n'), []]
        else: d[k][1].append(l.strip())
    return d
for src, dst, isprot in ((faa, out + '.proteins.faa', True), (fna, out + '.cds.fna', False)):
    d = fasta(src); changed = 0
    for tx, (cut5, nN) in trim5.items():
        if tx in d and nN:
            s = ''.join(d[tx][1]); n = nN // 3 if isprot else nN
            lead = s[:n]
            if (isprot and set(lead) <= {'X'}) or (not isprot and set(lead) <= {'N'}): d[tx][1] = [s[n:]]; changed += 1
            else: print(f'# {tx}: expected leading {"X" if isprot else "N"}x{n}, found {lead!r} — left unchanged', file=sys.stderr)
    with open(dst, 'w') as w:
        for k, (h, ss) in d.items():
            s = ''.join(ss); w.write(h + '\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
    print(f'{dst}: {changed} sequences trimmed', file=sys.stderr)
open(out + '.fixlog.tsv', 'w').write('\n'.join('\t'.join(r) for r in log) + '\n')
print(f'{len(log) - 1} transcripts flagged partial; 5-trimmed {len(trim5)}, 3-trimmed {len(trim3)}; CDS ids made unique', file=sys.stderr)
