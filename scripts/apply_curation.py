#!/usr/bin/env python3
"""Apply curation edits to a YaHS AGP and write the curated AGP + FASTA.
Edits (from join_qc.py / curation_contacts.py, see CURATION_WORKLIST_v3.md):
  breaks: split a scaffold at a component boundary (bg < 5 joins)
  joins : place scaffold X adjacent to an end of scaffold T; X's orientation is chosen from the Hi-C BAM so that
          the half of X with more contacts to T's terminal 1 Mb faces T.
Usage: apply_curation.py in.agp contigs.fa hic_scaffold.bam out_prefix"""
import sys, subprocess, collections
agp, ctgfa, bam, out = sys.argv[1:5]
BREAKS = [("scaffold_27", 5748000), ("scaffold_59", 927537), ("scaffold_119", 126000)]
JOINS = [("scaffold_52", "scaffold_25", "start"), ("scaffold_60", "scaffold_11", "start"), ("scaffold_44", "scaffold_7", "start")]
# ---- read AGP into ordered component lists (W rows only; gaps re-inserted on output) ----
scaf = collections.OrderedDict()
for l in open(agp):
    if l.startswith('#'): continue
    p = l.rstrip('\n').split('\t')
    if p[4] == 'W': scaf.setdefault(p[0], []).append([p[5], int(p[6]), int(p[7]), p[8], int(p[1]), int(p[2])])
L = {s: v[-1][5] for s, v in scaf.items()}
log = []
# ---- breaks ----
for s, pos in BREAKS:
    comps = scaf[s]; i = next(k for k, c in enumerate(comps) if c[5] == pos)
    left, right = comps[:i+1], comps[i+1:]
    scaf[s] = left; scaf[s + "_b"] = right
    log.append(f"BREAK {s} at {pos} ({comps[i][0]} | {comps[i+1][0]}): {s} keeps {len(left)} components, {s}_b gets {len(right)}")
# ---- orientation for joins from the BAM (contacts between X halves and T's terminal 1 Mb) ----
want = {x for x, _, _ in JOINS} | {t for _, t, _ in JOINS}
cnt = collections.Counter()
p = subprocess.Popen(["samtools", "view", "-q", "1", "-f", "64", "-F", "0x904", bam], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    f = line.split('\t', 9); c1, p1, c2, p2 = f[2], int(f[3]), f[6], int(f[7])
    if c2 == '=' or c1 == c2 or c1 not in want or c2 not in want: continue
    for X, T, end in JOINS:
        for (cx, px, ct, pt) in ((c1, p1, c2, p2), (c2, p2, c1, p1)):
            if cx == X and ct == T:
                inT = pt <= 1_000_000 if end == "start" else pt >= L[T] - 1_000_000
                if inT: cnt[(X, "first" if px < L[X] / 2 else "second")] += 1
for X, T, end in JOINS:
    a, b = cnt[(X, "first")], cnt[(X, "second")]
    # X is placed BEFORE T's start (so X's facing half must be its end) or AFTER T's end (facing half = its start)
    if end == "start":
        orient = "+" if b >= a else "-"          # '+' puts X's second half next to T's start
        newcomps = [[c[0], c[1], c[2], c[3]] for c in scaf[X]]
        if orient == "-": newcomps = [[c[0], c[1], c[2], "-" if c[3] == "+" else "+"] for c in reversed(newcomps)]
        scaf[T] = [c + [0, 0] for c in newcomps] + scaf[T]
    else:
        orient = "+" if a >= b else "-"          # '+' puts X's first half next to T's end
        newcomps = [[c[0], c[1], c[2], c[3]] for c in scaf[X]]
        if orient == "-": newcomps = [[c[0], c[1], c[2], "-" if c[3] == "+" else "+"] for c in reversed(newcomps)]
        scaf[T] = scaf[T] + [c + [0, 0] for c in newcomps]
    del scaf[X]
    log.append(f"JOIN {X} ({L[X]:,} bp) to {end} of {T}, orientation {orient} (contacts with T's terminal 1 Mb: first half {a}, second half {b})")
# ---- write AGP (100-bp U gaps between components) and FASTA ----
ctg = {}; name = None; seq = []
for l in open(ctgfa):
    if l.startswith('>'):
        if name: ctg[name] = ''.join(seq)
        name = l[1:].split()[0]; seq = []
    else: seq.append(l.strip())
ctg[name] = ''.join(seq)
comp_rc = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n'}
def rc(s): return ''.join(comp_rc.get(ch, 'N') for ch in reversed(s))
# order: by total length descending, renumber
def tlen(cs): return sum(c[2]-c[1]+1 for c in cs) + 100*(len(cs)-1)
order = sorted(scaf.items(), key=lambda kv: -tlen(kv[1]))
fa = open(out + ".fa", "w"); ag = open(out + ".agp", "w"); ag.write("##agp-version 2.1\n# curated from YaHS scaffolds; see CURATION_WORKLIST_v3.md\n")
namemap = open(out + ".names.tsv", "w"); namemap.write("curated_name\told_name\n")
for k, (old, cs) in enumerate(order, 1):
    new = f"scaffold_{k}"; namemap.write(f"{new}\t{old}\n"); pos = 1; part = 1; pieces = []
    for i, c in enumerate(cs):
        if i:
            ag.write(f"{new}\t{pos}\t{pos+99}\t{part}\tU\t100\tscaffold\tyes\tproximity_ligation\n"); pos += 100; part += 1; pieces.append('N'*100)
        ln = c[2]-c[1]+1; ag.write(f"{new}\t{pos}\t{pos+ln-1}\t{part}\tW\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\n"); pos += ln; part += 1
        s = ctg[c[0]][c[1]-1:c[2]]; pieces.append(rc(s) if c[3] == '-' else s)
    fa.write(f">{new}\n"); S = ''.join(pieces)
    for i in range(0, len(S), 100): fa.write(S[i:i+100] + "\n")
fa.close(); ag.close(); namemap.close()
open(out + ".edits.log", "w").write("\n".join(log) + "\n"); print("\n".join(log))
