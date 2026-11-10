#!/usr/bin/env python3
"""Curated mitogenome annotation of Mimp_v3_mito (18,133 bp), written as GenBank-style 5-column feature table (.tbl, for NGDC
GenBase) and as GFF3 (release / Data file 7). Starts from the MitoHiFi/MitoFinder calls (06_release/Mimp_v3_mito.gb) and applies
the corrections found on 2026-11-10 by checking every CDS against the sequence and the RefSeq annotation of the same species
(NC_011815.1), which also defines the conventions used (gene names trnX / 12S rRNA; incomplete stops via transl_except; ND3 frameshift):
  COX1  5422->5404: the true start is the in-frame GTG at 5404 (Val initiator read as Met; protein 516 aa, as in NC_011815.1);
        MitoFinder had started at the first ATG-less codon and the first GenBase upload therefore flagged it as 5' partial.
  ND4   10615->10255: true start ATG at 10255 (overlapping ND4L by 7 nt, as in the reference); 459 aa; ends with an incomplete
        T stop at 11632 completed by polyadenylation (transl_except), instead of running 5 nt into tRNA-His.
  ND3   join(9545..9718,9720..9894): the +1 frameshift at gene position 174 shared by turtles and birds (Mindell et al. 1998);
        116 aa, incomplete T stop at 9894 (transl_except).
  COX3  8693..9476 with incomplete T stop (transl_except), not 9475.
  ND1   2789..3759 with incomplete TA stop (transl_except), avoiding a 1-nt overlap with tRNA-Ile (reference convention).
  tRNA/rRNA get gene features (trnF … trnP; 12S rRNA / 16S rRNA) as requested by the GenBase curators; the four tRNA-Thr/tRNA-Pro
  units of the tandem array are annotated individually with a note; control region 16327..18133 as D-loop.
usage: 17b_mito_annotation_curated.py Mimp_v3_mito.fa out_prefix   -> out_prefix.tbl, out_prefix.gff3, out_prefix.proteins.faa"""
import sys
fa, out = sys.argv[1], sys.argv[2]
seqid = open(fa).readline()[1:].split()[0]; seq = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper(); L = len(seq)
CODE = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'M','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'W','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'*','AGG':'*','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
RC = str.maketrans('ACGT', 'TGCA')
POLYA = "TAA stop codon is completed by the addition of 3' A residues to the mRNA"
# (type, gene, product, strand, [(start,end),...], extra)   extra: te=(pos string) ; note ; start_codon
F = [
 ('tRNA','trnF','tRNA-Phe','+',[(1,70)],{}),
 ('rRNA','12S rRNA','12S ribosomal RNA','+',[(71,1041)],{}),
 ('tRNA','trnV','tRNA-Val','+',[(1042,1111)],{}),
 ('rRNA','16S rRNA','16S ribosomal RNA','+',[(1112,2710)],{}),
 ('tRNA','trnL2','tRNA-Leu','+',[(2711,2787)],{'note':'codons recognized: UUR'}),
 ('CDS','ND1','NADH dehydrogenase subunit 1','+',[(2789,3759)],{'te':'3758..3759','note':POLYA}),
 ('tRNA','trnI','tRNA-Ile','+',[(3760,3830)],{}),
 ('tRNA','trnQ','tRNA-Gln','-',[(3830,3900)],{}),
 ('tRNA','trnM','tRNA-Met','+',[(3900,3968)],{}),
 ('CDS','ND2','NADH dehydrogenase subunit 2','+',[(3969,5009)],{}),
 ('tRNA','trnW','tRNA-Trp','+',[(5024,5096)],{}),
 ('tRNA','trnA','tRNA-Ala','-',[(5098,5166)],{}),
 ('tRNA','trnN','tRNA-Asn','-',[(5168,5240)],{}),
 ('tRNA','trnC','tRNA-Cys','-',[(5266,5331)],{}),
 ('tRNA','trnY','tRNA-Tyr','-',[(5332,5402)],{}),
 ('CDS','COX1','cytochrome c oxidase subunit I','+',[(5404,6951)],{'note':'GTG start codon'}),
 ('tRNA','trnS2','tRNA-Ser','-',[(6943,7014)],{'note':'codons recognized: UCN'}),
 ('tRNA','trnD','tRNA-Asp','+',[(7015,7085)],{}),
 ('CDS','COX2','cytochrome c oxidase subunit II','+',[(7086,7772)],{}),
 ('tRNA','trnK','tRNA-Lys','+',[(7777,7849)],{}),
 ('CDS','ATP8','ATP synthase F0 subunit 8','+',[(7852,8019)],{}),
 ('CDS','ATP6','ATP synthase F0 subunit 6','+',[(8010,8693)],{}),
 ('CDS','COX3','cytochrome c oxidase subunit III','+',[(8693,9476)],{'te':'9476','note':POLYA}),
 ('tRNA','trnG','tRNA-Gly','+',[(9477,9544)],{}),
 ('CDS','ND3','NADH dehydrogenase subunit 3','+',[(9545,9718),(9720,9894)],{'te':'9894','note':'frameshift mechanism unknown (Mindell et al., 1998, Mol. Biol. Evol., 15:1568-1571); '+POLYA}),
 ('tRNA','trnR','tRNA-Arg','+',[(9895,9964)],{}),
 ('CDS','ND4L','NADH dehydrogenase subunit 4L','+',[(9965,10261)],{}),
 ('CDS','ND4','NADH dehydrogenase subunit 4','+',[(10255,11632)],{'te':'11632','note':POLYA}),
 ('tRNA','trnH','tRNA-His','+',[(11633,11701)],{}),
 ('tRNA','trnS1','tRNA-Ser','+',[(11702,11767)],{'note':'codons recognized: AGY'}),
 ('tRNA','trnL1','tRNA-Leu','+',[(11767,11837)],{'note':'codons recognized: CUN'}),
 ('CDS','ND5','NADH dehydrogenase subunit 5','+',[(11838,13619)],{}),
 ('CDS','ND6','NADH dehydrogenase subunit 6','-',[(13615,14139)],{}),
 ('tRNA','trnE','tRNA-Glu','-',[(14140,14207)],{}),
 ('CDS','CYTB','cytochrome b','+',[(14212,15348)],{}),
 ('tRNA','trnT','tRNA-Thr','+',[(15357,15426)],{'note':'tandem repeat unit 1 of 4 of the tRNA-Thr/tRNA-Pro array; the array is length-heteroplasmic (4 units in 44 %, 5 units in 52 % of HiFi reads)'}),
 ('tRNA','trnP','tRNA-Pro','-',[(15436,15505)],{'note':'tandem repeat unit 1 of 4'}),
 ('tRNA','trnT','tRNA-Thr','+',[(15613,15682)],{'note':'tandem repeat unit 2 of 4'}),
 ('tRNA','trnP','tRNA-Pro','-',[(15692,15761)],{'note':'tandem repeat unit 2 of 4'}),
 ('tRNA','trnT','tRNA-Thr','+',[(15922,15991)],{'note':'tandem repeat unit 3 of 4'}),
 ('tRNA','trnP','tRNA-Pro','-',[(16001,16070)],{'note':'tandem repeat unit 3 of 4'}),
 ('tRNA','trnT','tRNA-Thr','+',[(16178,16247)],{'note':'tandem repeat unit 4 of 4'}),
 ('tRNA','trnP','tRNA-Pro','-',[(16257,16326)],{'note':'tandem repeat unit 4 of 4'}),
 ('D-loop','','control region','+',[(16327,18133)],{'note':'control region'}),
]
# ---- checks: every CDS translates without internal stops; incomplete stops are T/TA
def cds_seq(parts, strand):
    s = ''.join(seq[a-1:b] for a, b in parts); return s.translate(RC)[::-1] if strand == '-' else s
def translate(s): return ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3))
prots = []
for typ, gene, prod, strand, parts, x in F:
    if typ != 'CDS': continue
    s = cds_seq(parts, strand); aa = translate(s); body = aa[:-1] if len(s) % 3 == 0 else aa
    assert '*' not in body, (gene, 'internal stop')
    if 'te' in x: assert len(s) % 3 in (1, 2) and s[len(s) - len(s) % 3:] in ('T', 'TA'), (gene, 'bad incomplete stop'); prot = aa
    else: assert len(s) % 3 == 0 and aa.endswith('*'), (gene, 'no stop'); prot = aa[:-1]
    if prot[0] != 'M': prot = 'M' + prot[1:]      # alternative initiator (GTG in COX1) is read as Met
    prots.append((gene, prot, len(s)))
# ---- .tbl (GenBank 5-column feature table)
def loc(a, b, strand): return (f'{a}\t{b}' if strand == '+' else f'{b}\t{a}')
tbl = [f'>Feature {seqid}']
for typ, gene, prod, strand, parts, x in F:
    key = typ if typ != 'D-loop' else 'D-loop'
    if typ in ('CDS', 'tRNA', 'rRNA'):
        tbl.append(loc(min(a for a, b in parts), max(b for a, b in parts), strand) + '\tgene'); tbl.append(f'\t\t\tgene\t{gene}')
    segs = parts if strand == '+' else parts[::-1]
    tbl.append(loc(*segs[0], strand) + f'\t{key}')
    for a, b in segs[1:]: tbl.append(loc(a, b, strand))
    if typ in ('CDS', 'tRNA', 'rRNA'): tbl.append(f'\t\t\tgene\t{gene}'); tbl.append(f'\t\t\tproduct\t{prod}')
    if typ == 'CDS':
        tbl.append('\t\t\tcodon_start\t1'); tbl.append('\t\t\ttransl_table\t2')
        if 'te' in x: tbl.append(f'\t\t\ttransl_except\t(pos:{x["te"]},aa:TERM)')
        if gene == 'ND3': tbl.append('\t\t\texception\tribosomal slippage')
    if 'note' in x and typ != 'D-loop': tbl.append(f'\t\t\tnote\t{x["note"]}')
    if typ == 'D-loop': tbl.append(f'\t\t\tnote\t{x["note"]}')
open(out + '.tbl', 'w').write('\n'.join(tbl) + '\n')
# ---- GFF3
g = ['##gff-version 3', f'##sequence-region {seqid} 1 {L}', f'{seqid}\tcurated\tregion\t1\t{L}\t.\t+\t.\tID={seqid};Name={seqid};Is_circular=true;genome=mitochondrion;mol_type=genomic DNA;organism=Manouria impressa']
n = 0
for typ, gene, prod, strand, parts, x in F:
    n += 1; gid = f'MIMP3_mt_g{n:02d}'; s, e = min(a for a, b in parts), max(b for a, b in parts)
    if typ == 'D-loop':
        g.append(f'{seqid}\tcurated\tD_loop\t{s}\t{e}\t.\t+\t.\tID=MIMP3_mt_dloop;Name=control region;Note=control region'); continue
    bt = {'CDS': 'protein_coding', 'tRNA': 'tRNA', 'rRNA': 'rRNA'}[typ]
    g.append(f'{seqid}\tcurated\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid};Name={gene};gene={gene};gene_biotype={bt}')
    cid = gid.replace('_g', {'CDS': '_cds', 'tRNA': '_t', 'rRNA': '_r'}[typ]); extra = ''
    if typ == 'CDS':
        extra = ';transl_table=2' + (f';transl_except=(pos:{x["te"]},aa:TERM)' if 'te' in x else '') + (';exception=ribosomal slippage' if gene == 'ND3' else '')
    note = f';Note={x["note"]}' if 'note' in x else ''
    run = 0
    for a, b in (parts if strand == '+' else parts[::-1]):
        phase = str((3 - run % 3) % 3) if typ == 'CDS' else '.'
        g.append(f'{seqid}\tcurated\t{typ}\t{a}\t{b}\t.\t{strand}\t{phase}\tID={cid};Parent={gid};Name={gene};product={prod};gene={gene}{extra}{note}'); run += b - a + 1
open(out + '.gff3', 'w').write('\n'.join(g) + '\n')
with open(out + '.proteins.faa', 'w') as w:
    for gene, prot, n in prots: w.write(f'>{seqid}_{gene} {len(prot)} aa; CDS {n} nt\n{prot}\n')
print('features:', len(F), '| CDS:', len(prots), '| proteins:', ', '.join(f'{g} {len(p)}' for g, p, n in prots))
