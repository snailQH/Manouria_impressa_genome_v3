#!/usr/bin/env python3
"""Unambiguous assembly statistics (the kind prior peer review asked for): sequences, total length, N50/N90,
longest, gap count and N bases; contig-level statistics obtained by splitting at runs of >=10 N.
Usage: asm_stats.py assembly.fa"""
import sys, re
L = []; name = None; seq = []
def flush():
    if name: L.append((name, ''.join(seq)))
for line in open(sys.argv[1]):
    if line.startswith('>'):
        flush(); name = line[1:].split()[0]; seq = []
    else: seq.append(line.strip())
flush()
def nstat(v, f):
    v = sorted(v, reverse=True); t = sum(v); a = 0
    for x in v:
        a += x
        if a >= t * f: return x
sc = [len(s) for _, s in L]; ctg = []; gaps = 0; nbases = 0
for _, s in L:
    parts = re.split(r'N{10,}', s.upper()); gaps += len(parts) - 1; nbases += s.upper().count('N'); ctg += [len(p) for p in parts if p]
print(f"sequences\t{len(sc)}\ntotal_bp\t{sum(sc)}\nscaffold_N50\t{nstat(sc,.5)}\nscaffold_N90\t{nstat(sc,.9)}\nlongest_scaffold\t{max(sc)}")
print(f"contigs(split_at_>=10N)\t{len(ctg)}\ncontig_N50\t{nstat(ctg,.5)}\ncontig_N90\t{nstat(ctg,.9)}\nlongest_contig\t{max(ctg)}\ngaps\t{gaps}\nN_bases\t{nbases}")
big = sorted(sc, reverse=True); print("top_sequences_Mb\t" + ",".join(f"{x/1e6:.1f}" for x in big[:40]))
