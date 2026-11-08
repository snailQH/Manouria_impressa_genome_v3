#!/usr/bin/env python3
"""Telomere check at chromosome ends: count vertebrate telomeric repeat (TTAGGG/CCCTAA) in the terminal windows.
usage: 09_telomere.py genome.fa [window_bp=50000] -> TSV on stdout (one row per chromosome end)."""
import sys, re, subprocess
fa = sys.argv[1]; W = int(sys.argv[2]) if len(sys.argv) > 2 else 50000
fai = [l.split("\t") for l in open(fa + ".fai")]
fwd, rev = re.compile("(?:TTAGGG){3,}"), re.compile("(?:CCCTAA){3,}")
print("sequence\tlength\tend\twindow_bp\ttelomere_repeat_bp\tlongest_array_bp\tstatus")
for name, L, *_ in fai:
    if not name.startswith("chr"): continue
    L = int(L)
    for end, reg in (("start", f"{name}:1-{min(W, L)}"), ("end", f"{name}:{max(1, L - W + 1)}-{L}")):
        s = "".join(subprocess.run(["samtools", "faidx", fa, reg], capture_output=True, text=True).stdout.split("\n")[1:]).upper()
        arr = [m.end() - m.start() for p in (fwd, rev) for m in p.finditer(s)]
        tot, mx = sum(arr), (max(arr) if arr else 0)
        print(f"{name}\t{L}\t{end}\t{W}\t{tot}\t{mx}\t{'telomere' if mx >= 300 else ('trace' if tot >= 60 else '-')}")
