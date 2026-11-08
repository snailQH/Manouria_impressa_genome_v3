#!/usr/bin/env python3
"""Name and orient the curated scaffolds for release.
  - chromosome-scale scaffolds (the largest N, N = expected haploid number) -> chr1..chrN by length
  - orientation: reverse-complement a chromosome if >50% of its aligned bases to the reference are on the minus strand
    (so that chromosomes read in the same direction as the related reference), else keep
  - all other scaffolds -> scaffold_<k> by length (unplaced)
Writes release FASTA + AGP (component = curated scaffold, so the curated AGP composes through) + a name/orientation table.
Usage: name_chromosomes.py curated.fa curated.agp aln.paf n_chromosomes out_prefix"""
import sys, collections
fa, agp, paf, N, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
seqs = collections.OrderedDict(); name = None; buf = []
for l in open(fa):
    if l.startswith('>'):
        if name: seqs[name] = ''.join(buf)
        name = l[1:].split()[0]; buf = []
    else: buf.append(l.strip())
seqs[name] = ''.join(buf)
order = sorted(seqs, key=lambda s: -len(seqs[s]))
# strand vote per scaffold from PAF (alignments >= 50 kb, mapq >= 30), weighted by aligned length
plus = collections.Counter(); minus = collections.Counter(); best = collections.defaultdict(collections.Counter)
for l in open(paf):
    p = l.split('\t')
    if int(p[10]) < 50000 or int(p[11]) < 30: continue
    (plus if p[4] == '+' else minus)[p[0]] += int(p[10]); best[p[0]][p[5]] += int(p[10])
comp = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n'}
def rc(s): return ''.join(comp.get(c, 'N') for c in reversed(s))
tab = open(out + '.names.tsv', 'w'); tab.write("release_name\tcurated_name\tlength\torientation\tbest_reference_chr\taligned_bp\n")
fo = open(out + '.fa', 'w'); ag = open(out + '.agp', 'w'); ag.write("##agp-version 2.1\n")
for i, s in enumerate(order, 1):
    L = len(seqs[s]); new = f"chr{i}" if i <= N else f"scaffold_{i - N}"
    flip = minus[s] > plus[s]
    seq = rc(seqs[s]) if flip else seqs[s]
    ref = best[s].most_common(1)[0] if best[s] else ('-', 0)
    tab.write(f"{new}\t{s}\t{L}\t{'-' if flip else '+'}\t{ref[0]}\t{ref[1]}\n")
    ag.write(f"{new}\t1\t{L}\t1\tW\t{s}\t1\t{L}\t{'-' if flip else '+'}\n")
    fo.write(f">{new}\n"); [fo.write(seq[j:j+100] + "\n") for j in range(0, L, 100)]
print(f"wrote {out}.fa/.agp/.names.tsv: {N} chromosomes + {len(order)-N} unplaced; flipped: {sum(1 for s in order[:N] if minus[s] > plus[s])} chromosomes")
