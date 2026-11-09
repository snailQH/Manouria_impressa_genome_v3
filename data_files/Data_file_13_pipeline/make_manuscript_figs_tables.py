#!/usr/bin/env python3
"""Builds the manuscript figures and tables that can be produced from the release files.
Inputs: docs/genome_paper_new_submission/results/inputs/ (collected by collect_manuscript_inputs.sh).
Outputs: docs/genome_paper_new_submission/figures/ (PNG, 300 dpi) and supplementary/ + manuscript/tables/ (TSV/MD).
Annotation-dependent panels (repeat landscape, gene summary) are added by the same script once their inputs exist."""
import os, re, csv, collections, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUB = os.path.join(BASE, 'docs', 'genome_paper_new_submission'); IN = os.path.join(SUB, 'results', 'inputs')
FIGROOT = os.path.join(SUB, 'figures'); FIG = os.path.join(FIGROOT, 'panels'); BMC = os.path.join(FIGROOT, 'bmc'); LONG = os.path.join(FIGROOT, 'longform')
SUP = os.path.join(SUB, 'supplementary'); TAB = os.path.join(SUB, 'manuscript', 'tables')
for d in (FIG, BMC, LONG, SUP, TAB): os.makedirs(d, exist_ok=True)
# panels/: single-panel source images (regenerated here). bmc/ and longform/: the composites each manuscript version actually uses.
# A panel appears in only one place per manuscript version; composites are built from panels/ so nothing is stored twice by hand.
P = lambda f: os.path.join(IN, f)
INK, INK2, HIFI, HIC, ANN, GREY = '#1c1f1a', '#5e655b', '#2b5d8c', '#8a4b1f', '#3f6b3a', '#c7c9bf'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})

def fai(f): return [(r[0], int(r[1])) for r in csv.reader(open(P(f)), delimiter='\t')]
def nx(lengths, x):
    L = sorted(lengths, reverse=True); t = sum(L) * x / 100; c = 0
    for l in L:
        c += l
        if c >= t: return l
def tsv(f): return list(csv.DictReader(open(P(f)), delimiter='\t'))
def fmt(n): return f'{n:,}'

asm = {'Mimp_v3': fai('Mimp_v3.fa.fai'), 'hap1': fai('Mimp_v3_hap1.fa.fai'), 'hap2': fai('Mimp_v3_hap2.fa.fai'), 'v2.4': fai('v2.4.fa.fai')}
G = sum(l for _, l in asm['Mimp_v3'])

# ---------- Figure 2D: Nx curves and cumulative length ----------
col = {'Mimp_v3': HIFI, 'hap1': ANN, 'hap2': '#95c48d', 'v2.4': GREY}
for variant, keys, out in (('main', ('Mimp_v3', 'hap1', 'hap2'), 'nx_cumulative.png'),):
  fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
  for k in keys:
    v = asm[k]; k = 'preliminary .1' if k == 'v2.4' else k; colk = col['v2.4'] if k == 'preliminary .1' else col[k]
    L = sorted((l for _, l in v), reverse=True); xs = list(range(1, 101)); a1.plot(xs, [nx(L, x) / 1e6 for x in xs], color=colk, lw=1.6, label=k)
    cum, s = [], 0
    for l in L[:60]: s += l; cum.append(s / 1e9)
    a2.step(range(1, len(cum) + 1), cum, where='post', color=colk, lw=1.6, label=k)
  a1.set_xlabel('x (%)'); a1.set_ylabel('Nx scaffold length (Mb)'); a1.axvline(50, color=GREY, lw=0.8, ls=':'); a1.legend(frameon=False, fontsize=8)
  a2.set_xlabel('Scaffolds, longest first'); a2.set_ylabel('Cumulative length (Gb)'); a2.axhline(G / 1e9, color=GREY, lw=0.8, ls=':'); a2.set_ylim(0, 2.6)
  fig.tight_layout(); fig.savefig(os.path.join(FIG, out), dpi=300); plt.close(fig)

# ---------- chromosome table S1 + Figure 3A' chromosome summary ----------
names = {r['release_name']: r for r in tsv('Mimp_v3.names.tsv')}
comp = {r['sequence']: r for r in tsv('agp_components.tsv')}
tel = collections.defaultdict(dict)
for r in tsv('Mimp_v3.telomeres.tsv'): tel[r['sequence']][r['end']] = r['status']
v24 = {r['v3_chr']: r['v2.4_sequences_mapping_>=50%'] for r in tsv('v24_names_per_v3_chr.tsv')}
snv = {h: {r['sequence']: float(r['SNV_per_kb']) for r in tsv(f'snv_density_{h}_vs_primary.tsv')} for h in ('hap1', 'hap2')}
chrs = [(n, l) for n, l in asm['Mimp_v3'] if n.startswith('chr')]
with open(os.path.join(SUP, 'TableS1_chromosomes.tsv'), 'w') as f:
    f.write('chromosome\tlength_bp\tcontigs\tgaps\torientation_vs_GopEvg\tG_evgoodei_homolog\taligned_bp_to_homolog\ttelomere_start\ttelomere_end\tSNV_per_kb_hap1\tSNV_per_kb_hap2\tv2.4_sequences\n')
    for n, l in chrs:
        r = names[n]; c = comp.get(n, {'contigs': '1', 'gaps': '0'})
        f.write('\t'.join([n, str(l), c['contigs'], c['gaps'], r['orientation'], r['best_reference_chr'], r['aligned_bp'], tel[n].get('start', '-'), tel[n].get('end', '-'),
                           f"{snv['hap1'].get(n, 0):.2f}", f"{snv['hap2'].get(n, 0):.2f}", v24.get(n, '-')]) + '\n')
# ideogram-style figure
cagp = [r for r in csv.reader(open(P('curated.agp')), delimiter='\t') if not r[0].startswith('#')]
fig, ax = plt.subplots(figsize=(7.2, 5.2)); y = 0
for n, l in chrs:
    ax.add_patch(Rectangle((0, y - 0.32), l / 1e6, 0.64, facecolor='#dde8f2', edgecolor=HIFI, lw=0.6))
    cur = names[n]['curated_name']
    for r in cagp:
        if r[0] == cur and r[4] in ('N', 'U'):
            x = int(r[1]) / 1e6; x = l / 1e6 - x if names[n]['orientation'] == '-' else x
            ax.plot([x] * 2, [y - 0.32, y + 0.32], color=HIC, lw=0.9)
    for end, x in (('start', 0), ('end', l / 1e6)):
        st = tel[n].get(end, '-')
        if st == 'telomere': ax.plot(x, y, marker='o', ms=4.5, color=ANN, zorder=5)
        elif st == 'trace': ax.plot(x, y, marker='o', ms=3, mfc='none', mec=ANN, zorder=5)
    frac = int(names[n]['aligned_bp']) / l if names[n]['aligned_bp'].isdigit() else 0
    lab = names[n]['best_reference_chr'].replace('GopEvg_', 'GopEvg ') + ('' if frac >= 0.5 else f' ({100*frac:.0f} %)') if frac >= 0.1 else 'largely unaligned (repeat-rich)'
    ax.text(l / 1e6 + 4, y, lab, va='center', fontsize=6.5, color=INK2)
    y += 1
ax.set_yticks(range(len(chrs))); ax.set_yticklabels([n for n, _ in chrs], fontsize=7); ax.invert_yaxis(); ax.set_xlabel('Position (Mb)'); ax.set_xlim(0, 420)
ax.plot([], [], color=HIC, lw=1, label='scaffold gap (contig junction)'); ax.plot([], [], 'o', ms=4.5, color=ANN, label='telomeric array ≥ 300 bp'); ax.plot([], [], 'o', ms=3, mfc='none', mec=ANN, label='telomeric trace')
ax.legend(frameon=False, fontsize=7, loc='lower right'); ax.set_title('Mimp_v3 chromosomes: contig junctions, telomeres, G. evgoodei homolog', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'chromosome_summary.png'), dpi=300); plt.close(fig)

# ---------- panel: per-chromosome Merqury QV ----------
if os.path.exists(P('merqury_Mimp_v3_per_sequence.qv')):
    qvc = {l.split()[0]: float(l.split()[3]) for l in open(P('merqury_Mimp_v3_per_sequence.qv')) if l.split()[0].startswith('chr')}
    fig, ax = plt.subplots(figsize=(7.2, 2.6)); xs = range(len(chrs))
    ax.bar(list(xs), [qvc.get(n, 0) for n, _ in chrs], color=HIFI); ax.axhline(float(open(P('merqury_Mimp_v3.qv')).read().split()[3]), color=GREY, lw=0.8, ls=':')
    ax.set_xticks(list(xs)); ax.set_xticklabels([n.replace('chr', '') for n, _ in chrs], fontsize=7); ax.set_xlabel('Chromosome'); ax.set_ylabel('Merqury QV'); ax.set_ylim(50, 92)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'per_chromosome_qv.png'), dpi=300); plt.close(fig)

# ---------- panel: haplotype divergence per chromosome ----------
fig, ax = plt.subplots(figsize=(7.2, 2.8)); xs = range(len(chrs)); w = 0.4
ax.bar([x - w / 2 for x in xs], [snv['hap1'].get(n, 0) for n, _ in chrs], w, color=ANN, label='hap1 vs primary')
ax.bar([x + w / 2 for x in xs], [snv['hap2'].get(n, 0) for n, _ in chrs], w, color='#95c48d', label='hap2 vs primary')
ax.set_xticks(list(xs)); ax.set_xticklabels([n.replace('chr', '') for n, _ in chrs], fontsize=7); ax.set_xlabel('Chromosome'); ax.set_ylabel('SNVs per kb'); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'haplotype_divergence_per_chromosome.png'), dpi=300); plt.close(fig)
with open(os.path.join(SUP, 'TableS8_haplotype_divergence.tsv'), 'w') as f:
    f.write(open(P('hap_divergence_summary.tsv')).read()); f.write('\n# per-chromosome SNV density (hap vs primary): see TableS1_chromosomes.tsv\n')

# ---------- Supplementary figure: mitogenome repeat-array copy number ----------
h = [(int(r[0]), int(r[1])) for r in csv.reader(open(P('mito_array_copynumber_hist.tsv')), delimiter=' ')]
fig, ax = plt.subplots(figsize=(4.2, 2.6)); tot = sum(c for _, c in h)
ax.bar([b for b, _ in h], [100 * c / tot for _, c in h], width=90, color=HIC); ax.set_xlim(-700, 1300)
ax.set_xlabel('Net indel within the tRNA-Thr/Pro array per spanning read (bp)'); ax.set_ylabel('Reads (%)')
ax.text(0, 46, 'four units\n(consensus)', ha='center', fontsize=7, color=INK2); ax.text(300, 52, 'five units', ha='center', fontsize=7, color=INK2); ax.set_ylim(0, 60)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'mito_array_heteroplasmy.png'), dpi=300); plt.close(fig)

# ---------- long-form Figure 2 (A GenomeScope, B Nx/cumulative, C per-chromosome QV, D synteny dot plot) ----------
# Complementary to the Data Note's Fig 1 (Hi-C map, chromosome summary, repeat landscape, Merqury spectrum): no panel is shared.
import matplotlib.image as mpimg
def composite(panels, out, size=(7.2, 6.4), upper=True):
    fig = plt.figure(figsize=size); gs = fig.add_gridspec(2, 2, hspace=0.12, wspace=0.08)
    for i, (lab, f) in enumerate(panels):
        ax = fig.add_subplot(gs[i // 2, i % 2]); ax.axis('off'); ax.text(0.0, 1.02, lab, transform=ax.transAxes, fontsize=12, fontweight='bold')
        if f and os.path.exists(f): ax.imshow(mpimg.imread(f))
        else: ax.text(0.5, 0.5, '[[PENDING]]', ha='center', va='center', fontsize=9, color=INK2)
    fig.savefig(out, dpi=300, bbox_inches='tight'); plt.close(fig)

# ---------- Table 2 ----------
def qc_compleasm(f):
    t = open(P(f)).read(); g = lambda k: re.search(k + r':([\d.]+)%', t); return f"S {g('S').group(1)} D {g('D').group(1)} F {g('F').group(1)} M {g('M').group(1)}" if g('S') else '—'
def qc_busco(f):
    m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)', open(P(f)).read()); return f'C {m.group(1)} [S {m.group(2)}, D {m.group(3)}] F {m.group(4)} M {m.group(5)}' if m else '—'
qv = open(P('merqury_Mimp_v3.qv')).read().split()[3]; cpl = open(P('merqury_Mimp_v3.completeness.stats')).read().split()[-1]
rows = [['Statistic', 'Mimp_v3 (reference)', 'hap1', 'hap2', 'preliminary version .1']]
def stat(name, fn): rows.append([name] + [fn(k) for k in ('Mimp_v3', 'hap1', 'hap2', 'v2.4')])
stat('Sequences', lambda k: fmt(len(asm[k]))); stat('Total length (bp)', lambda k: fmt(sum(l for _, l in asm[k])))
stat('Scaffold N50 (bp)', lambda k: fmt(nx([l for _, l in asm[k]], 50))); stat('Scaffold N90 (bp)', lambda k: fmt(nx([l for _, l in asm[k]], 90)))
stat('Longest sequence (bp)', lambda k: fmt(max(l for _, l in asm[k])))
def stats(f):
    d = dict(l.rstrip('\n').split('\t')[:2] for l in open(P(f)) if '\t' in l); return (fmt(int(d['contigs(split_at_>=10N)'])), fmt(int(d['contig_N50'])), fmt(int(d['gaps'])))
contigs = {'Mimp_v3': stats('Mimp_v3.stats.tsv'), 'hap1': stats('Mimp_v3_hap1.stats.tsv'), 'hap2': stats('Mimp_v3_hap2.stats.tsv'), 'v2.4': ('502', '86,871,185', '100')}
rows.append(['Contigs'] + [contigs[k][0] for k in ('Mimp_v3', 'hap1', 'hap2', 'v2.4')]); rows.append(['Contig N50 (bp)'] + [contigs[k][1] for k in ('Mimp_v3', 'hap1', 'hap2', 'v2.4')]); rows.append(['Gaps'] + [contigs[k][2] for k in ('Mimp_v3', 'hap1', 'hap2', 'v2.4')])
rows.append(['Chromosome-scale sequences', f'26 ({100*sum(l for n,l in chrs)/G:.1f} % of length)', '—', '—', '32 (26 chromosomes + 6 fragments)'])
hq = {l.split()[0]: l.split()[3] for l in open(P('merqury_haplotypes.qv'))}; hc = {l.split()[0]: l.split()[4] for l in open(P('merqury_haplotypes.completeness.stats'))}
rows.append(['Merqury QV / k-mer completeness (%)', f'{float(qv):.2f} / {float(cpl):.2f}', f"{float(hq['asm1']):.2f} / {float(hc['asm1']):.2f}", f"{float(hq['asm2']):.2f} / {float(hc['asm2']):.2f} (both haplotypes: QV {float(hq['Both']):.2f}, completeness {float(hc['both']):.2f})", '67.43 / 95.43'])
rows.append(['compleasm vertebrata_odb12 (n = 3,390) (%)', qc_compleasm('compleasm_Mimp_v3_vertebrata_odb12.txt'), qc_compleasm('compleasm_hap1_scaffolds_vertebrata_odb12.txt'), qc_compleasm('compleasm_hap2_scaffolds_vertebrata_odb12.txt'), qc_compleasm('compleasm_v2.4_previous_vertebrata_odb12.txt')])
rows.append(['BUSCO 5.8.0 vertebrata_odb10 (n = 3,354) (%)', qc_busco('busco_Mimp_v3.txt'), qc_busco('busco_hap1_scaffolds.txt'), qc_busco('busco_hap2_scaffolds.txt'), qc_busco('busco_v2.4_previous.txt')])
def ebp(k, q):
    L=[l for _, l in asm[k]]; cn50 = {'Mimp_v3': 86871185, 'hap1': 36595271, 'hap2': 28880847, 'v2.4': 86871185}[k]
    return f"{int(math.log10(cn50))}.{int(math.log10(nx(L, 50)))}.Q{int(float(q))}"
rows.append(['EBP quality notation', ebp('Mimp_v3', qv), ebp('hap1', hq['asm1']), ebp('hap2', hq['asm2']), ebp('v2.4', 67.43)])
rows.append(['Accession', 'GWHHIVU00000000.2', 'GWHKKMH00000000.1', 'GWHKKMI00000000.1', 'GWHHIVU00000000.1'])
seq = [['Library', 'Platform', 'Reads / pairs', 'Bases (Gb)', 'Mean / N50 read length (bp)', 'Coverage', 'Accession'],
       ['PacBio HiFi (liver)', 'PacBio Revio, 4 SMRT cells, ccs 7.0.0', '3,618,287 reads (3,618,282 after adapter screen)', '68.55', '18,945 / 18,784', '28.2× of assembly; 30× of k-mer estimate', 'CRR1746698 (CRA024571)'],
       ['Hi-C (liver, MboI)', 'MGI DNBSEQ, 2 × 150 bp', '32,770,909 pairs', '9.83', '150', '20,760,724 valid pairs (63.4 %)', 'CRR1746699 (CRA024571)'],
       ['RNA-seq heart (IT3)', 'MGI DNBSEQ-G400 [[confirm]], 2 × 150 bp', '22,319,011 pairs', '6.70', '150', '—', 'CRR1747397 (CRA024586)'],
       ['RNA-seq lung (IT4)', 'MGI DNBSEQ-G400 [[confirm]], 2 × 150 bp', '22,448,804 pairs', '6.73', '150', '—', 'CRR1747398 (CRA024586)']]
def write_table(name, blocks):
    with open(os.path.join(TAB, name + '.md'), 'w') as f:
        for title, rws in blocks:
            f.write(f'**{title}**\n\n| ' + ' | '.join(rws[0]) + ' |\n|' + '---|' * len(rws[0]) + '\n')
            for r in rws[1:]: f.write('| ' + ' | '.join(r) + ' |\n')
            f.write('\n')
    with open(os.path.join(TAB, name + '.tsv'), 'w') as f:
        for title, rws in blocks:
            f.write(f'# {title}\n'); [f.write('\t'.join(r) + '\n') for r in rws]
main_rows = [r[:4] for r in rows]
write_table('Table2_sequencing_assembly', [('Table 2a. Sequencing data (counted from the raw files)', seq), ('Table 2b. Assembly statistics', main_rows)])
with open(os.path.join(SUP, 'TableS4b_preliminary_version_statistics.tsv'), 'w') as f:
    f.write('# Assembly statistics of the preliminary Genome Warehouse version .1 next to Mimp_v3 (same reads; contigs identical)\n')
    for r in rows: f.write('\t'.join([r[0], r[1], r[4]]) + '\n')

# ---------- Table 1 (software) ----------
sw = tsv('software_versions.tsv'); step_names = {'1': 'read preparation', '2': 'contig assembly', '3': 'Hi-C mapping/scaffolding', '3b': 'haplotype scaffolding', '3c': 'contact maps', '3e': 'synteny', '3f': 'Hi-C statistics', '4': 'assembly QC', '4c': 'coverage', '5': 'repeat library', '11': 'contamination / mitogenome', '12': 'repeat masking', '13': 'gene prediction', '14': 'functional / ncRNA annotation', '15': 'release files'}
t1 = [['Tool', 'Version', 'Used for', 'Non-default parameters']]
params = {'hifiasm': '--h1/--h2 (Hi-C mode)', 'chromap': '--preset hic -q 1 --remove-pcr-duplicates', 'yahs': '-e GATC', 'HapHiC': 'pipeline --RE GATC, 26 groups', 'minimap2': '-x asm20 (synteny), map-hifi (coverage), asm5 (haplotypes, v2.4)', 'meryl': 'k = 21', 'genomescope2': '-k 21 -p 2', 'compleasm': 'vertebrata_odb12, sauropsida_odb12', 'busco': 'vertebrata_odb10, offline', 'RepeatModeler': '-LTRStruct', 'RepeatMasker': '-xsmall; species library then -species tetrapoda', 'BRAKER3': '--bam (HISAT2) --prot_seq OrthoDB v12 Vertebrata --busco_lineage vertebrata_odb10', 'eggnog-mapper': '-m diamond --tax_scope 7742 --go_evidence non-electronic', 'MitoHiFi': '-f/-g NC_011815.1 -o 2', 'FCS-GX': '--tax-id 74924', 'FCS-adaptor': '--euk', 'samtools': '', 'seqkit': '', 'pretextmap': '', 'pretextsnapshot': '', 'juicer_tools': 'pre', 'merqury': 'k = 21', 'assembly-stats': '', 'diamond': '', 'bedtools': '', 'gfatools': ''}
for r in sw:
    steps = ', '.join(step_names.get(x.strip(), x.strip()) for x in r['step'].split(','))
    t1.append([r['tool'], r['version'], steps, params.get(r['tool'], '')])
t1 += [['tRNAscan-SE', '2.0.12', 'tRNA annotation', '-E -H --detail; EukHighConfidenceFilter'], ['Barrnap', '0.9', 'rRNA annotation', '--kingdom euk'], ['paftools.js', 'minimap2 2.28', 'haplotype divergence', 'call -l 10000 -L 50000'], ['TSEBRA', 'BRAKER 3.0.8', 'gene-set merge (compared, not used for the final set)', 'default'], ['bedtools', '2.31.1', 'RNA-seq support per gene; soft-mask union', 'multicov; maskfasta -soft']]
write_table('Table1_software', [('Table 1. Software, versions and non-default parameters', t1)])

# ---------- Supplementary tables copied / derived ----------
import shutil
shutil.copy(P('Mimp_v3.hic_stats.tsv'), os.path.join(SUP, 'TableS3_hic_statistics.tsv'))
shutil.copy(P('v2.4_to_v3_correspondence.tsv'), os.path.join(SUP, 'TableS4_v2.4_to_v3_correspondence.tsv'))
shutil.copy(P('Mimp_v3.telomeres.tsv'), os.path.join(SUP, 'TableS5_telomeres.tsv'))
with open(os.path.join(SUP, 'TableS6_haplotype_to_primary.tsv'), 'w') as f:
    for h in ('hap1', 'hap2'):
        f.write(f'# {h}\n'); f.write(open(P(f'Mimp_v3_{h}_to_primary.tsv')).read())
# S2 curation edits: edits log + after-curation scores of the new junctions
with open(os.path.join(SUP, 'TableS2_curation_edits.tsv'), 'w') as f:
    f.write('# Edits applied by apply_curation.py (v3_curated.edits.log)\n'); f.write(open(P('curation_edits.log')).read())
    f.write('\n# Join scores after curation (join_qc_v3_curated.tsv: ratio = observed/expected for contiguous sequence; bg = observed/inter-scaffold background)\n'); f.write(open(P('join_qc_after_curation.tsv')).read())
# S7 coverage anomalies ≥ 50 kb
if os.path.exists(P('fcs_gx_actions_applied.tsv')): shutil.copy(P('fcs_gx_actions_applied.tsv'), os.path.join(SUP, 'TableS11_fcs_gx_removed_sequences.tsv'))
with open(os.path.join(SUP, 'TableS7_coverage_zero_runs.tsv'), 'w') as f:
    f.write('sequence\tstart\tend\tlength_bp\tinterpretation\n')
    for r in csv.reader(open(P('coverage_zero_runs.bed')), delimiter='\t'):
        L = int(r[2]) - int(r[1])
        if L >= 50000:
            note = 'heterozygous tandem duplication (reads present at MAPQ 0)' if r[0] in ('chr1', 'scaffold_28') else 'zero coverage at MAPQ ≥ 1; see text'
            f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{L}\t{note}\n')
# ---------- repeats (when masking has run) ----------
if os.path.exists(P('repeat_landscape_Mimp_v3.png')):
    shutil.copy(P('repeat_landscape_Mimp_v3.png'), os.path.join(FIG, 'repeat_landscape.png'))
    t = open(P('repeatmasker_custom.tbl')).read(); rows = [['Repeat class', 'Elements', 'Length (bp)', 'Genome (%)']]
    for line in t.splitlines():
        m = re.match(r'\s*([A-Za-z0-9/\-]+(?: [A-Za-z\-]+)*):?\s+(\d+)\s+(\d+) bp\s+([\d.]+) %', line)
        if m and m.group(1) not in ('CRE/SLACS',): rows.append([m.group(1), fmt(int(m.group(2))), fmt(int(m.group(3))), m.group(4)])
    m2 = re.search(r'bases masked:\s+(\d+) bp \(\s*([\d.]+) %\)', open(P('repeatmasker_dfam.tbl')).read())
    rows.append(['Additional masking, Dfam 3.8 Tetrapoda pass', '', fmt(int(m2.group(1))), m2.group(2)])
    with open(os.path.join(SUP, 'TableS9_repeat_classes.tsv'), 'w') as f:
        f.write('# RepeatMasker 4.1.7-p1 with the RepeatModeler2 library (genome.fa.tbl); second pass with Dfam Tetrapoda on the masked output\n')
        for r in rows: f.write('\t'.join(r) + '\n')
# ---------- Supplementary Fig S4: haplotype Hi-C maps ----------
if os.path.exists(P('hic_map_Mimp_v3_hap1.png')) and os.path.exists(P('hic_map_Mimp_v3_hap2.png')):
    fig = plt.figure(figsize=(7.2, 3.8)); gs = fig.add_gridspec(1, 2, wspace=0.05)
    for i, h in enumerate(('hap1', 'hap2')):
        ax = fig.add_subplot(gs[0, i]); ax.imshow(mpimg.imread(P(f'hic_map_Mimp_v3_{h}.png'))); ax.axis('off'); ax.text(0.0, 1.02, 'AB'[i], transform=ax.transAxes, fontsize=12, fontweight='bold'); ax.set_title(f'Mimp_v3_{h}', fontsize=9, loc='right')
    fig.savefig(os.path.join(FIG, 'haplotype_hic_maps.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
# ---------- Supplementary Table S12: synteny vs G. evgoodei ----------
if os.path.exists(P('Mimp_v3_vs_GopEvg_blocks.tsv')):
    with open(os.path.join(SUP, 'TableS12_synteny_GopEvg.tsv'), 'w') as f:
        for name, title in (('composition', 'A. Reference-chromosome composition of each Mimp_v3 chromosome (synteny blocks >= 1 Mb)'), ('breakpoints', 'B. Candidate structural differences and Hi-C / contig evidence on the Mimp_v3 side'), ('blocks', 'C. Synteny blocks')):
            f.write(f'# {title}\n'); f.write(open(P(f'Mimp_v3_vs_GopEvg_{name}.tsv')).read()); f.write('\n')
# ---------- BMC Genomic Data Fig 1 composite ----------
hic = P('hic_map_Mimp_v3_pub.png') if os.path.exists(P('hic_map_Mimp_v3_pub.png')) else P('hic_map_Mimp_v3.png')
if os.path.exists(os.path.join(FIG, 'chromosome_summary.png')):
    fig = plt.figure(figsize=(7.2, 7.6)); gs = fig.add_gridspec(2, 2, hspace=0.08, wspace=0.06)
    panels = [('a', hic), ('b', os.path.join(FIG, 'chromosome_summary.png')), ('c', os.path.join(FIG, 'repeat_landscape.png')), ('d', P('merqury_Mimp_v3_spectra-cn.png'))]
    for i, (lab, f) in enumerate(panels):
        ax = fig.add_subplot(gs[i // 2, i % 2]); ax.imshow(mpimg.imread(f)); ax.axis('off'); ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=13, fontweight='bold')
    fig.savefig(os.path.join(BMC, 'Fig1_assembly_overview.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
composite([('A', P('genomescope_linear_plot.png')), ('B', os.path.join(FIG, 'nx_cumulative.png')), ('C', os.path.join(FIG, 'per_chromosome_qv.png')), ('D', P('dotplot_Mimp_v3_vs_GopEvg.png'))], os.path.join(LONG, 'Fig1_assembly_validation.png'))
# long-form Fig 3 composite: A chromosome summary, B repeat landscape, D mitogenome (C = gene summary, added when the gene set is final)
fig = plt.figure(figsize=(7.2, 7.2)); gs = fig.add_gridspec(2, 2, hspace=0.08, wspace=0.06)
pan3 = [('A', os.path.join(FIG, 'gene_summary.png') if os.path.exists(os.path.join(FIG, 'gene_summary.png')) else None), ('B', os.path.join(FIG, 'haplotype_divergence_per_chromosome.png')), ('C', P('mito_annotation.png')), ('D', os.path.join(FIG, 'mito_array_heteroplasmy.png'))]
for i, (lab, f) in enumerate(pan3):
    ax = fig.add_subplot(gs[i // 2, i % 2]); ax.axis('off'); ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=13, fontweight='bold')
    if f and os.path.exists(f): ax.imshow(mpimg.imread(f))
    else: ax.text(0.5, 0.5, '[[PENDING: gene annotation summary]]', ha='center', va='center', fontsize=9, color=INK2)
fig.savefig(os.path.join(LONG, 'Fig2_annotation_haplotypes_mitogenome.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
# long-form supplementary figures = copies of panels (so the longform/ folder is complete on its own)
for f in os.listdir(LONG):
    if f.startswith('FigS'): os.remove(os.path.join(LONG, f))
for src, dst in (('haplotype_hic_maps.png', 'FigS1_haplotype_hic_maps.png'),):
    sp = os.path.join(FIG, src) if not src.startswith('..') else P(os.path.basename(src))
    if os.path.exists(sp): shutil.copy(sp, os.path.join(LONG, dst))
# ---------- gene annotation: Table 3, Table S10, gene-summary panel ----------
if os.path.exists(P('gene_set_stats.tsv')):
    gs_ = dict(l.rstrip('\n').split('\t')[:2] for l in open(P('gene_set_stats.tsv')) if '\t' in l)
    bp = qc_busco('busco_protein_final.txt') if os.path.exists(P('busco_protein_final.txt')) else '—'
    eg = open(P('eggnog_rate.txt')).read().split() if os.path.exists(P('eggnog_rate.txt')) else None
    rs = open(P('rnaseq_support_final.txt')).read().strip() if os.path.exists(P('rnaseq_support_final.txt')) else ''
    m_rs = re.search(r'>=10 reads (\d+)/(\d+) \(([\d.]+) %\); >=1 read (\d+) \(([\d.]+) %\)', rs)
    t3 = [['Feature', 'Value'],
          ['Protein-coding genes', fmt(int(gs_['genes']))], ['Transcripts', fmt(int(gs_['transcripts']))], ['Transcripts per gene', gs_['mean_transcripts_per_gene']],
          ['Mean gene length (bp)', fmt(int(gs_['mean_gene_length_bp']))], ['Mean CDS length (bp)', fmt(int(gs_['mean_CDS_length_bp']))], ['Mean exons per transcript', gs_['mean_exons_per_transcript']],
          ['Mean exon / intron length (bp)', f"{gs_['mean_exon_length_bp']} / {fmt(int(gs_['mean_intron_length_bp']))}"], ['Median protein length (aa)', gs_['median_protein_length_aa']], ['Single-exon transcripts', fmt(int(gs_['single_exon_transcripts']))],
          ['BUSCO vertebrata_odb10, protein mode (n = 3,354) (%)', bp],
          ['Genes with eggNOG-mapper annotation (%)', f"{fmt(int(eg[0]))} ({100*int(eg[0])/int(eg[1]):.1f}); GO {fmt(int(eg[2]))}, KEGG {fmt(int(eg[3]))}, Pfam {fmt(int(eg[4]))}" if eg else '—'],
          ['Genes with RNA-seq support (heart + lung)', f"{fmt(int(m_rs.group(1)))} ({m_rs.group(3)} %) with ≥ 10 reads; {fmt(int(m_rs.group(4)))} ({m_rs.group(5)} %) with ≥ 1 read" if m_rs else '—'],
          ['tRNA genes (tRNAscan-SE high-confidence / all predictions)', '3,309 / 24,407'], ['rRNA genes or fragments (Barrnap)', '18S 92; 5.8S 79; 28S 99; 5S 181'],
          ['Repeats (soft-masked, union of two passes) (%)', '48.60']]
    write_table('Table3_annotation', [('Table 3. Annotation statistics (final evidence-supported gene set)', t3)])
    shutil.copy(P('gene_set_stats.tsv'), os.path.join(SUP, 'TableS10_gene_set_statistics.tsv'))
    if os.path.exists(P('busco_protein_candidates.txt')):
        with open(os.path.join(SUP, 'TableS10_gene_set_statistics.tsv'), 'a') as f: f.write('\n# BUSCO-protein of candidate gene sets (decision: highest C with D < 4 %, then evidence filter)\n' + open(P('busco_protein_candidates.txt')).read())
    # gene-summary panel: genes per chromosome (bar) + gene density
    if os.path.exists(P('gene_density.tsv')):
        gd = {r[0]: (int(r[1]), float(r[3])) for r in csv.reader(open(P('gene_density.tsv')), delimiter='\t') if r and r[0] != 'sequence' and r[1].isdigit()}
        fig, ax = plt.subplots(figsize=(7.2, 2.8)); xs = range(len(chrs))
        ax.bar(list(xs), [gd.get(n, (0, 0))[0] for n, _ in chrs], color=ANN); ax.set_ylabel('Protein-coding genes'); ax.set_xlabel('Chromosome')
        ax.set_xticks(list(xs)); ax.set_xticklabels([n.replace('chr', '') for n, _ in chrs], fontsize=7)
        ax2 = ax.twinx(); ax2.plot(list(xs), [gd.get(n, (0, 0))[1] for n, _ in chrs], 'o-', color=HIC, ms=3, lw=1); ax2.set_ylabel('Genes per Mb', color=HIC); ax2.spines['right'].set_visible(True)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, 'gene_summary.png'), dpi=300); plt.close(fig)
if os.path.exists(os.path.join(FIG, 'gene_summary.png')):
    composite([('A', os.path.join(FIG, 'gene_summary.png')), ('B', os.path.join(FIG, 'haplotype_divergence_per_chromosome.png')), ('C', P('mito_annotation.png')), ('D', os.path.join(FIG, 'mito_array_heteroplasmy.png'))], os.path.join(LONG, 'Fig2_annotation_haplotypes_mitogenome.png'), size=(7.2, 7.2))
print('panels:', sorted(os.listdir(FIG))); print('bmc:', sorted(os.listdir(BMC))); print('longform:', sorted(os.listdir(LONG))); print('supplementary:', sorted(os.listdir(SUP))); print('tables:', sorted(os.listdir(TAB)))
