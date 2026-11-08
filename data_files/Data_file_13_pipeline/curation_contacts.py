#!/usr/bin/env python3
"""Quantitative Hi-C checks for the curation worklist, from the scaffold-level Hi-C BAM.
Reads read-1 records (flag 0x40, MAPQ>=1) and uses RNEXT/PNEXT for the mate, so each pair is counted once.
  (1) for each (small scaffold X, candidate chromosome T): pairs X-T binned along T (5-Mb bins) and along X,
      plus X's pairs with every other chromosome-scale scaffold, per Mb, as background;
  (2) for a scaffold with a suspected internal seam: contact between the two halves vs within-half at the
      same genomic distance, and vs. that scaffold's contact with other chromosomes.
Usage: curation_contacts.py hic.bam scaffolds.fa.fai out_prefix"""
import sys, subprocess, collections
bam, fai, out = sys.argv[1:4]
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
CHR = [f"scaffold_{i}" for i in range(1, 27)]
import os
CASES = eval(os.environ["CASES"]) if "CASES" in os.environ else {"scaffold_22": "scaffold_3", "scaffold_44": "scaffold_7", "scaffold_52": "scaffold_25", "scaffold_60": "scaffold_11",
         "scaffold_35": "scaffold_1", "scaffold_39": "scaffold_18", "scaffold_58": "scaffold_9", "scaffold_17": None, "scaffold_24": None, "scaffold_27": None}
SEAM = eval(os.environ["SEAM"]) if "SEAM" in os.environ else {"scaffold_12": 0.58}
B = 5_000_000
pairs_XT = collections.defaultdict(lambda: collections.Counter())   # (X,T) -> bin on T
pairs_Xchr = collections.defaultdict(collections.Counter)            # X -> chr -> n
seam = {s: collections.Counter() for s in SEAM}                      # intra counts by (binA,binB) 1-Mb
inter_chr = collections.Counter()                                    # pairs between different chromosome-scale scaffolds (per chr) for background
p = subprocess.Popen(["samtools", "view", "-q", "1", "-f", "64", "-F", "0x904", bam], stdout=subprocess.PIPE, text=True)
n = 0
for line in p.stdout:
    f = line.split("\t", 9); c1, p1, c2, p2 = f[2], int(f[3]), f[6], int(f[7])
    if c2 == "=": c2 = c1
    n += 1
    for X, Y in ((c1, c2), (c2, c1)):
        if X in CASES and Y != X:
            pairs_Xchr[X][Y] += 1
            T = CASES[X]
            if T and Y == T:
                posT = p2 if X == c1 else p1; posX = p1 if X == c1 else p2
                pairs_XT[(X, T)][posT // B] += 1; pairs_XT[(X, "selfpos")][posX // 500_000] += 1
    if c1 in SEAM and c2 == c1:
        a, b = sorted((p1 // 1_000_000, p2 // 1_000_000)); seam[c1][(a, b)] += 1
    if c1 != c2 and c1 in CHR and c2 in CHR:
        inter_chr[c1] += 1; inter_chr[c2] += 1
o = open(out + "_contacts.txt", "w")
def w(s): print(s); o.write(s + "\n")
w(f"read pairs scanned: {n:,}")
tot_chr_len = sum(L[c] for c in CHR)
for X, T in CASES.items():
    lx = L[X]; byc = pairs_Xchr[X]; tot_other = sum(v for c, v in byc.items() if c in CHR)
    w(f"\n== {X} ({lx/1e6:.2f} Mb): pairs with chromosome-scale scaffolds = {tot_other:,}  ({tot_other/(lx/1e6):.0f} per Mb of X)")
    top = byc.most_common(5)
    w("   top partners (pairs, pairs per Mb of partner, enrichment vs partner's share of genome):")
    for c, v in top:
        if c in L: w(f"     {c:12s} {v:7,d}  {v/(L[c]/1e6):7.1f}/Mb   enrich {v/tot_other/(L[c]/tot_chr_len):5.2f}x" if tot_other else "")
    if T:
        bins = pairs_XT[(X, T)]; nb = L[T] // B + 1; vals = [bins.get(i, 0) for i in range(nb)]
        sp = pairs_XT[(X, "selfpos")]; nx = L[X] // 500_000 + 1; sv = [sp.get(i, 0) for i in range(nx)]
        w(f"   where on {X} (500-kb bins) the {T} contacts come from: {sv}")
        w(f"   pairs with {T} along {T} in 5-Mb bins (first 4 | ... | last 4): {vals[:4]} ... {vals[-4:]}   max bin {max(vals)} at {vals.index(max(vals))*5} Mb; median bin {sorted(vals)[len(vals)//2]}")
for s, frac in SEAM.items():
    cut = int(L[s] * frac / 1_000_000); M = seam[s]; nb = L[s] // 1_000_000 + 1
    def mean_at(d, side):  # mean pairs between 1-Mb bins at distance d, restricted to side: 'A' (both<cut), 'B' (both>=cut), 'X' (across)
        vals = []
        for a in range(nb - d):
            b = a + d
            if side == 'A' and b < cut: vals.append(M.get((a, b), 0))
            elif side == 'B' and a >= cut: vals.append(M.get((a, b), 0))
            elif side == 'X' and a < cut <= b: vals.append(M.get((a, b), 0))
        return sum(vals) / len(vals) if vals else float('nan')
    w(f"\n== {s} seam test at {cut} Mb (fraction {frac}): mean pairs per 1-Mb bin pair at distance d")
    w("   d(Mb)   within-A   within-B   across-seam   across/within")
    for d in (2, 4, 6, 8, 10, 15, 20):
        A, Bv, X = mean_at(d, 'A'), mean_at(d, 'B'), mean_at(d, 'X'); wi = (A + Bv) / 2
        w(f"   {d:5d}   {A:8.1f}   {Bv:8.1f}   {X:11.1f}   {X/wi if wi else float('nan'):6.2f}")
    xs = sum(v for (a, b), v in M.items() if a < cut <= b); inter = inter_chr[s]
    w(f"   total pairs across the seam: {xs:,}; {s}'s pairs with OTHER chromosomes: {inter:,} (over {tot_chr_len/1e6:.0f} Mb) -> background {inter/(tot_chr_len/1e6):.1f} per Mb of partner; across-seam per Mb of part B: {xs/((L[s]-cut*1e6)/1e6):.1f}")
o.close()
