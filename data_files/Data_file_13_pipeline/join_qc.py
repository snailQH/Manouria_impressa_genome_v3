#!/usr/bin/env python3
"""Score every contig join in a scaffold assembly from Hi-C read pairs.

For each join (boundary between consecutive AGP components inside a scaffold) the observed number of read
pairs linking the W kb immediately before the boundary to the W kb immediately after it is compared with the
number expected for sequence that is truly contiguous. Expectation comes from the genome-wide distance-decay
curve measured INSIDE hifiasm contigs (pairs per 100-kb bin pair as a function of separation), so no join
contributes to its own expectation. ratio = observed / expected:
   ~1        contiguous (correct join)
   << 1      the two sides do not interact -> misjoin (or an unmappable/centromeric stretch: check coverage)
A second statistic, 'bg', is the ratio observed / inter-scaffold background for the same two windows, so a
misjoin reads as bg ~ 1 and a correct join as bg >> 1.
Also scores user-specified candidate joins between scaffold ends (for curation decisions and for verifying
applied edits).
Usage: join_qc.py hic.bam scaffolds.agp scaffolds.fa.fai out.tsv [candidates.tsv]
   candidates.tsv: scaffoldA  endA(start|end)  scaffoldB  endB(start|end)   (one per line)
"""
import sys, subprocess, collections, math
bam, agp, fai, out = sys.argv[1:5]; cand = sys.argv[5] if len(sys.argv) > 5 else None
BIN = 100_000; W = 20          # 100-kb bins; window = 20 bins = 2 Mb each side
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
# contigs and joins from AGP
comp = collections.defaultdict(list)
for l in open(agp):
    if l.startswith('#'): continue
    p = l.rstrip('\n').split('\t')
    if p[4] == 'W': comp[p[0]].append((int(p[1]), int(p[2]), p[5], p[8]))
joins = []
for s, cs in comp.items():
    for i in range(1, len(cs)):
        joins.append((s, cs[i-1][1], cs[i-1][2], cs[i][2]))        # (scaffold, boundary pos, left contig, right contig)
# contig interval lookup per scaffold (for "inside a contig" test)
def contig_at(s, pos):
    for a, b, n, o in comp[s]:
        if a <= pos <= b: return n
    return None
# --- pass over read pairs: per-scaffold sparse 100-kb matrices (cis) + inter-scaffold totals ---
cis = collections.defaultdict(collections.Counter); trans = collections.Counter(); trans_pairs = 0
p = subprocess.Popen(["samtools", "view", "-q", "1", "-f", "64", "-F", "0x904", bam], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    f = line.split('\t', 9); c1, p1, c2, p2 = f[2], int(f[3]), f[6], int(f[7])
    if c2 == '=': c2 = c1
    if c1 == c2:
        a, b = sorted((p1 // BIN, p2 // BIN)); cis[c1][(a, b)] += 1
    else:
        trans[c1] += 1; trans[c2] += 1; trans_pairs += 1
# --- distance decay inside contigs (bins fully inside one contig, both ends in the same contig) ---
decay = collections.Counter(); ndec = collections.Counter()
for s, M in cis.items():
    if s not in comp: continue
    for (a, b), v in M.items():
        d = b - a
        if d == 0 or d > 2 * W: continue
        ca = contig_at(s, a * BIN + BIN // 2); cb = contig_at(s, b * BIN + BIN // 2)
        if ca and ca == cb: decay[d] += v; ndec[d] += 1
exp_at = {d: decay[d] / ndec[d] for d in decay if ndec[d] > 50}
genome = sum(L.values()); bg_per_binpair = trans_pairs * 2 / ((genome / BIN) ** 2)   # inter-scaffold pairs per bin pair (rough, symmetric)
def score(sA, posA, sideA, sB, posB, sideB):
    """observed pairs between the W bins on sideA of posA in sA and the W bins on sideB of posB in sB,
    expected = sum over bin pairs of decay(dist if contiguous)."""
    binsA = [posA // BIN - 1 - i for i in range(W)] if sideA == 'left' else [posA // BIN + i for i in range(W)]
    binsB = [posB // BIN + i for i in range(W)] if sideB == 'right' else [posB // BIN - 1 - i for i in range(W)]
    binsA = [b for b in binsA if 0 <= b <= L[sA] // BIN]; binsB = [b for b in binsB if 0 <= b <= L[sB] // BIN]
    obs = 0; exp = 0.0; n = 0
    for i, a in enumerate(binsA):
        for j, b in enumerate(binsB):
            d = i + j + 1
            if sA == sB: obs += cis[sA].get((min(a, b), max(a, b)), 0)
            exp += exp_at.get(d, exp_at[max(exp_at)]); n += 1
    return obs, exp, n
# cross-scaffold observed counts need a second pass only for candidates; do it cheaply by re-reading pairs for those scaffolds
cands = []
if cand:
    for l in open(cand):
        if l.strip() and not l.startswith('#'):
            a, ea, b, eb = l.split()[:4]; cands.append((a, ea, b, eb))
cross = collections.defaultdict(collections.Counter)
if cands:
    want = {c[0] for c in cands} | {c[2] for c in cands}
    p = subprocess.Popen(["samtools", "view", "-q", "1", "-f", "64", "-F", "0x904", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t', 9); c1, p1, c2, p2 = f[2], int(f[3]), f[6], int(f[7])
        if c2 == '=' or c1 == c2: continue
        if c1 in want and c2 in want:
            cross[(c1, c2)][(p1 // BIN, p2 // BIN)] += 1; cross[(c2, c1)][(p2 // BIN, p1 // BIN)] += 1
o = open(out, 'w'); o.write("scaffold\tboundary_bp\tleft_contig\tright_contig\tobserved\texpected\tratio_obs_exp\tbg_ratio\tverdict\n")
res = []
for s, pos, lc, rc in joins:
    obs, exp, n = score(s, pos, 'left', s, pos, 'right'); bg = bg_per_binpair * n
    r = obs / exp if exp else float('nan'); b = obs / bg if bg else float('nan')
    v = 'OK' if r >= 0.5 else ('LOW' if r >= 0.2 else 'BREAK?')
    res.append((s, pos, lc, rc, obs, exp, r, b, v)); o.write(f"{s}\t{pos}\t{lc}\t{rc}\t{obs}\t{exp:.0f}\t{r:.2f}\t{b:.1f}\t{v}\n")
for a, ea, b, eb in cands:
    posA = L[a] if ea == 'end' else 0; posB = L[b] if eb == 'end' else 0
    sideA = 'left' if ea == 'end' else 'right'; sideB = 'right' if eb == 'start' else 'left'
    binsA = [posA // BIN - 1 - i for i in range(W)] if sideA == 'left' else [i for i in range(W)]
    binsB = [i for i in range(W)] if sideB == 'right' else [posB // BIN - 1 - i for i in range(W)]
    binsA = [x for x in binsA if 0 <= x <= L[a] // BIN]; binsB = [x for x in binsB if 0 <= x <= L[b] // BIN]
    obs = sum(cross[(a, b)].get((x, y), 0) for x in binsA for y in binsB)
    exp = sum(exp_at.get(i + j + 1, exp_at[max(exp_at)]) for i in range(len(binsA)) for j in range(len(binsB))); n = len(binsA) * len(binsB)
    bg = bg_per_binpair * n; r = obs / exp if exp else float('nan'); bb = obs / bg if bg else float('nan')
    v = 'JOIN-SUPPORTED' if r >= 0.5 else ('WEAK' if r >= 0.2 else 'NOT-SUPPORTED')
    o.write(f"CANDIDATE:{a}:{ea}\t-\t{a}\t{b}:{eb}\t{obs}\t{exp:.0f}\t{r:.2f}\t{bb:.1f}\t{v}\n"); res.append((f"CAND {a}:{ea}+{b}:{eb}", 0, '', '', obs, exp, r, bb, v))
o.close()
rs = sorted([x[6] for x in res if not str(x[0]).startswith('CAND') and not math.isnan(x[6])])
print(f"joins scored: {len(rs)}; ratio median {rs[len(rs)//2]:.2f}, 10th pct {rs[len(rs)//10]:.2f}, min {rs[0]:.2f}; decay bins used: {len(exp_at)}")
print("flagged (ratio < 0.5):")
for x in res:
    if (str(x[0]).startswith('CAND')) or x[6] < 0.5: print(f"  {x[0]:40s} {x[1]:>11} {x[2]:>12} {x[3]:>12}  obs {x[4]:6d}  exp {x[5]:7.0f}  ratio {x[6]:5.2f}  bg {x[7]:6.1f}  {x[8]}")
