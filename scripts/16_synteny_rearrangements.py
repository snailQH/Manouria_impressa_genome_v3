#!/usr/bin/env python3
"""Synteny blocks between Mimp_v3 chromosomes and G. evgoodei, candidate structural differences, and the Hi-C
evidence for each candidate breakpoint in OUR assembly.

A difference to the reference is only reported as a feature of the M. impressa assembly when the Mimp_v3 side is
physically supported: (i) the breakpoint lies inside a single HiFi contig, and/or (ii) Hi-C contact across the
breakpoint is at the level expected for contiguous sequence (ratio ~1; join_qc statistic) — i.e. the assembly is
continuous there. Direction of change (fusion in M. impressa vs fission in G. evgoodei etc.) is NOT inferred.

usage: 16_synteny_rearrangements.py aln.paf hic.bam release.fa.fai names.tsv curated.agp out_prefix [min_block_bp=1000000]
"""
import sys, subprocess, collections
paf, bam, fai, names_tsv, cagp, outp = sys.argv[1:7]; MINB = int(sys.argv[7]) if len(sys.argv) > 7 else 1_000_000
BIN = 100_000; W = 20
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
chrs = [c for c in L if c.startswith('chr')]

# ---------- synteny blocks ----------
aln = collections.defaultdict(list)
for l in open(paf):
    p = l.split('\t')
    if int(p[10]) < 50_000 or int(p[11]) < 30 or p[0] not in L: continue
    aln[p[0]].append((int(p[2]), int(p[3]), p[5], p[4], int(p[7]), int(p[8])))
blocks = {}
for c, A in aln.items():
    A.sort(); B = []
    for qs, qe, t, st, ts, te in A:
        if B and B[-1][2] == t and B[-1][3] == st and qs - B[-1][1] < 3_000_000 and (abs(ts - B[-1][5]) < 3_000_000 if st == '+' else abs(B[-1][4] - te) < 3_000_000):
            b = B[-1]; B[-1] = (b[0], max(b[1], qe), t, st, min(b[4], ts), max(b[5], te), b[6] + (qe - qs))
        else: B.append((qs, qe, t, st, ts, te, qe - qs))
    blocks[c] = [b for b in B if b[1] - b[0] >= MINB]
with open(outp + '_blocks.tsv', 'w') as f:
    f.write('chromosome\tq_start\tq_end\tblock_bp\tref_chr\tstrand\tr_start\tr_end\taligned_bp\n')
    for c in chrs:
        for b in blocks.get(c, []): f.write(f'{c}\t{b[0]}\t{b[1]}\t{b[1]-b[0]}\t{b[2]}\t{b[3]}\t{b[4]}\t{b[5]}\t{b[6]}\n')

# ---------- candidate breakpoints between adjacent blocks ----------
cands = []   # (chr, pos, type, description)
for c in chrs:
    B = blocks.get(c, [])
    for i in range(1, len(B)):
        a, b = B[i-1], B[i]; pos = (a[1] + b[0]) // 2
        if a[2] != b[2]: cands.append((c, pos, 'inter-chromosomal', f'{a[2]} ({(a[1]-a[0])/1e6:.1f} Mb) | {b[2]} ({(b[1]-b[0])/1e6:.1f} Mb)'))
        elif a[3] != b[3]: cands.append((c, pos, 'strand switch (inversion candidate)', f'{a[2]} {a[3]}→{b[3]}'))
        else:
            gap = (b[4] - a[5]) if a[3] == '+' else (a[4] - b[5])
            if abs(gap) > 5_000_000: cands.append((c, pos, 'order jump (translocation candidate)', f'{a[2]} ref gap {gap/1e6:.1f} Mb'))
# chromosome composition
comp = collections.defaultdict(collections.Counter)
for c in chrs:
    for b in blocks.get(c, []): comp[c][b[2]] += b[6]

# ---------- contig map in release coordinates (curated AGP + orientation) ----------
rel = {}; orient = {}
for l in open(names_tsv):
    p = l.rstrip('\n').split('\t')
    if p[0] != 'release_name': rel[p[1]] = p[0]; orient[p[0]] = p[3]
contigs = collections.defaultdict(list)
for l in open(cagp):
    if l.startswith('#'): continue
    p = l.rstrip('\n').split('\t')
    if p[4] == 'W' and p[0] in rel:
        c = rel[p[0]]; s, e = int(p[1]), int(p[2])
        if orient.get(c) == '-': s, e = L[c] - e + 1, L[c] - s + 1
        contigs[c].append((s, e, p[5]))
def contig_at(c, pos):
    for s, e, n in contigs[c]:
        if s <= pos <= e: return n
    return None

# ---------- Hi-C: cis matrices + within-contig decay (as join_qc.py) ----------
cis = collections.defaultdict(collections.Counter); trans_pairs = 0
pr = subprocess.Popen(['samtools', 'view', '-q', '1', '-f', '64', '-F', '0x904', bam], stdout=subprocess.PIPE, text=True)
for line in pr.stdout:
    f = line.split('\t', 9); c1, p1, c2, p2 = f[2], int(f[3]), f[6], int(f[7])
    if c2 == '=': c2 = c1
    if c1 == c2:
        a, b = sorted((p1 // BIN, p2 // BIN)); cis[c1][(a, b)] += 1
    else: trans_pairs += 1
decay = collections.Counter(); ndec = collections.Counter()
for c in chrs:
    for (a, b), v in cis[c].items():
        d = b - a
        if d == 0 or d > 2 * W: continue
        ca, cb = contig_at(c, a * BIN + BIN // 2), contig_at(c, b * BIN + BIN // 2)
        if ca and ca == cb: decay[d] += v; ndec[d] += 1
exp_at = {d: decay[d] / ndec[d] for d in decay if ndec[d] > 50}
genome = sum(L.values()); bg = trans_pairs * 2 / ((genome / BIN) ** 2)
def across(c, pos):
    k = pos // BIN; obs = exp = 0.0
    for i in range(max(0, k - W), k):
        for j in range(k, min(L[c] // BIN + 1, k + W)):
            obs += cis[c].get((i, j), 0); exp += exp_at.get(j - i, 0)
    return obs, exp, (obs / exp if exp else float('nan')), (obs / (bg * W * W) if bg else float('nan'))

with open(outp + '_breakpoints.tsv', 'w') as f:
    f.write('chromosome\tposition\ttype\tdetail\tinside_contig\thic_pairs_across_2Mb\texpected_if_contiguous\tratio\tbg_ratio\tverdict\n')
    seen = set()
    for c, pos, typ, det in cands:
        if (c, pos) in seen: continue
        seen.add((c, pos)); ctg = contig_at(c, pos); obs, exp, r, b = across(c, pos)
        # Physical evidence for the Mimp_v3 side: a single HiFi contig spanning the breakpoint is decisive; at a contig
        # junction the join must pass the curation thresholds (ratio >= 0.5 or bg >= 11). A low ratio inside a contig
        # (centromere-like contact decay) does not contradict contiguity.
        if ctg: verdict = 'continuous in Mimp_v3: one HiFi contig spans the breakpoint; Hi-C across = %.2f of contiguous expectation, %.0fx background' % (r, b)
        elif r >= 0.5 or b >= 11: verdict = 'continuous in Mimp_v3: contig junction supported by Hi-C (ratio %.2f, %.0fx background)' % (r, b)
        else: verdict = 'NOT supported at a contig junction (ratio %.2f, %.0fx background) — possible misjoin, do not report' % (r, b)
        f.write(f'{c}\t{pos}\t{typ}\t{det}\t{ctg or "junction"}\t{obs:.0f}\t{exp:.0f}\t{r:.2f}\t{b:.1f}\t{verdict}\n')
with open(outp + '_composition.tsv', 'w') as f:
    f.write('chromosome\tlength\tref_chromosomes_(aligned_Mb_in_blocks_>=%dMb)\tn_blocks\tstrand_switches\tinter_chromosomal_junctions\n' % (MINB // 1_000_000))
    for c in chrs:
        B = blocks.get(c, []); sw = sum(1 for i in range(1, len(B)) if B[i][2] == B[i-1][2] and B[i][3] != B[i-1][3]); ic = sum(1 for i in range(1, len(B)) if B[i][2] != B[i-1][2])
        f.write(f"{c}\t{L[c]}\t{'; '.join(f'{k.replace(chr(95), chr(32))} {v/1e6:.1f}' for k, v in comp[c].most_common())}\t{len(B)}\t{sw}\t{ic}\n")
print('blocks:', sum(len(v) for v in blocks.values()), 'candidates:', len(cands), 'decay bins:', len(exp_at))
