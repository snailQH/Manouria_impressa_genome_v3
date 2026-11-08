#!/usr/bin/env python3
"""Repeat landscape from RepeatMasker .out files (divergence column = % substitutions from consensus, uncorrected).
usage: 12_repeat_landscape.py genome_size_bp out_png table_tsv file1.out [file2.out ...] [--tbl file.tbl]
With --tbl, each class histogram is scaled so its total equals the RepeatMasker .tbl class total (non-overlapping, authoritative).
Bins divergence in 1 % steps per class (LTR, LINE, SINE, DNA, RC, Unknown, Other); masked bp per bin as % of genome."""
import sys, collections
args = sys.argv[1:]; tbl = None
if '--tbl' in args: i = args.index('--tbl'); tbl = args[i+1]; args = args[:i] + args[i+2:]
G = int(args[0]); png, tsv, files = args[1], args[2], args[3:]
cls = ['LTR', 'LINE', 'SINE', 'DNA', 'RC', 'Unknown', 'Other']
land = collections.defaultdict(lambda: [0.0]*51); total = collections.Counter(); seen = set()
for f in files:
    for l in open(f):
        p = l.split()
        if len(p) < 15 or not p[0].isdigit(): continue
        div = float(p[1]); s, e = int(p[5]), int(p[6]); c = p[10].split('/')[0]
        if c in ('Simple_repeat', 'Low_complexity', 'Satellite', 'rRNA', 'tRNA', 'snRNA', 'scRNA', 'srpRNA'): continue
        key = (p[4], s, e)
        if key in seen: continue
        seen.add(key); c = c if c in cls else 'Other'
        land[c][min(int(div), 50)] += e - s + 1; total[c] += e - s + 1
if tbl:
    import re
    t = open(tbl).read(); key = {'SINE': r'SINEs:\s+\d+\s+(\d+) bp', 'LINE': r'LINEs:\s+\d+\s+(\d+) bp', 'LTR': r'LTR elements:\s+\d+\s+(\d+) bp', 'DNA': r'DNA transposons\s+\d+\s+(\d+) bp', 'RC': r'Rolling-circles\s+\d+\s+(\d+) bp', 'Unknown': r'Unclassified:\s+\d+\s+(\d+) bp'}
    for c, k in key.items():
        m = re.search(k, t)
        if m and total[c]: f = int(m.group(1)) / total[c]; land[c] = [v * f for v in land[c]]; total[c] = int(m.group(1))
with open(tsv, 'w') as o:
    o.write('class\t' + '\t'.join(f'div_{i}' for i in range(51)) + '\ttotal_bp\tpct_genome\n')
    for c in cls: o.write(c + '\t' + '\t'.join(f'{100*v/G:.4f}' for v in land[c]) + f'\t{total[c]}\t{100*total[c]/G:.2f}\n')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
col = {'LTR': '#8a4b1f', 'LINE': '#2b5d8c', 'SINE': '#6fa8dc', 'DNA': '#3f6b3a', 'RC': '#9a6200', 'Unknown': '#9aa09a', 'Other': '#d0d0d0'}
fig, ax = plt.subplots(figsize=(8, 4)); bottom = [0.0]*51; x = list(range(51))
for c in cls:
    y = [100*v/G for v in land[c]]; ax.bar(x, y, bottom=bottom, color=col[c], width=1, label=f'{c} ({100*total[c]/G:.1f} %)'); bottom = [b+v for b, v in zip(bottom, y)]
ax.set_xlabel('Divergence from consensus (%)'); ax.set_ylabel('Genome (%)'); ax.set_xlim(-0.5, 50.5); ax.legend(frameon=False, fontsize=8)
ax.set_title('Repeat landscape, Mimp_v3'); fig.tight_layout(); fig.savefig(png, dpi=200)
print('interspersed repeats:', f'{100*sum(total.values())/G:.2f} % of genome;', ', '.join(f'{c} {100*total[c]/G:.2f} %' for c in cls))
