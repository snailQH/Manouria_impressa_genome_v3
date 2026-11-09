#!/usr/bin/env python3
"""Data Note DOCX (BMC Genomic Data format with `bmc`; long form otherwise) from the Markdown single source of truth.
  docs/genome_paper_new_submission/manuscript/MANUSCRIPT_DATA_NOTE.md  (text; [n] citations)
  docs/genome_paper_new_submission/manuscript/references_verified.md   (numbered, verified references)
  docs/genome_paper_new_submission/manuscript/tables/*.md              (Tables 1-2)
  docs/genome_paper_new_submission/figures/*.png                       (figures)
Output: docs/genome_paper_new_submission/manuscript/Manouria_impressa_genome_v3_DataNote.docx
Placeholders [[...]] are kept verbatim and highlighted yellow so nothing pending can be mistaken for a result."""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from refnumber import Renumberer
from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX, WD_BREAK
from docx.shared import Inches, Pt

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUB = os.path.join(BASE, 'docs', 'genome_paper_new_submission'); MAN = os.path.join(SUB, 'manuscript'); FIG = os.path.join(SUB, 'figures', 'longform'); BMCFIG = os.path.join(SUB, 'figures', 'bmc'); RES = os.path.join(SUB, 'results', 'inputs')
OUT = os.path.join(MAN, 'Manouria_impressa_genome_v3_DataNote.docx')
INLINE = re.compile(r'(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+?`|\[\[.+?\]\])')

def styles(doc):
    n = doc.styles['Normal']; n.font.name = 'Times New Roman'; n.font.size = Pt(12); n.paragraph_format.line_spacing = 2.0; n.paragraph_format.space_after = Pt(0)
    t = doc.styles.add_style('ManuscriptTitle', WD_STYLE_TYPE.PARAGRAPH); t.font.name = 'Times New Roman'; t.font.size = Pt(16); t.font.bold = True; t.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER; t.paragraph_format.space_after = Pt(12)
    for name, size in (('Heading 1', 14), ('Heading 2', 13), ('Heading 3', 12)):
        h = doc.styles[name]; h.font.name = 'Times New Roman'; h.font.size = Pt(size); h.font.bold = True; h.font.color.rgb = None; h.paragraph_format.space_before = Pt(12); h.paragraph_format.space_after = Pt(6)
    c = doc.styles.add_style('Caption2', WD_STYLE_TYPE.PARAGRAPH); c.font.name = 'Times New Roman'; c.font.size = Pt(10); c.paragraph_format.line_spacing = 1.15; c.paragraph_format.space_after = Pt(8)

def rich(par, text, size=None):
    for chunk in INLINE.split(text):
        if not chunk: continue
        if chunk.startswith('[[') and chunk.endswith(']]'):
            r = par.add_run(chunk); r.font.highlight_color = WD_COLOR_INDEX.YELLOW; r.bold = True
        elif chunk.startswith('**') and chunk.endswith('**') and len(chunk) > 4: r = par.add_run(chunk[2:-2]); r.bold = True
        elif chunk.startswith('*') and chunk.endswith('*') and len(chunk) > 2: r = par.add_run(chunk[1:-1]); r.italic = True
        elif chunk.startswith('`') and chunk.endswith('`'): r = par.add_run(chunk[1:-1]); r.font.name = 'Courier New'
        else: r = par.add_run(chunk)
        if size: r.font.size = Pt(size)

def md_table(doc, lines, font=9):
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in lines if not re.match(r'^\s*\|[\s:|-]+\|\s*$', l)]
    ncol = max(len(r) for r in rows); t = doc.add_table(rows=0, cols=ncol); t.style = 'Table Grid'
    for i, cells in enumerate(rows):
        w = t.add_row().cells
        for j, c in enumerate(cells + [''] * (ncol - len(cells))):
            w[j].text = ''; p = w[j].paragraphs[0]; p.paragraph_format.line_spacing = 1.0; rich(p, c, font)
            if i == 0:
                for r in p.runs: r.bold = True
    doc.add_paragraph()

def references():
    """Numbered entries from references_verified.md, verification tags stripped."""
    refs = {}
    for l in open(os.path.join(MAN, 'references_verified.md'), encoding='utf-8'):
        m = re.match(r'^(\d+)\.\s+(.*)', l)
        if m: refs[int(m.group(1))] = m.group(2).strip()
    return refs

def main(src='MANUSCRIPT_DATA_NOTE.md', out=OUT, figs=None, author_placeholder=True, table_subs=None, extra_tables=True):
    doc = Document(); styles(doc)
    for s in doc.sections: s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Inches(1)
    md = open(os.path.join(MAN, src), encoding='utf-8').read().split('\n')
    # title = first H1; skip the draft-conventions block (lines before the first '---')
    title = md[0].lstrip('# ').strip(); doc.add_paragraph(title, style='ManuscriptTitle')
    if author_placeholder: p = doc.add_paragraph(); rich(p, '[[AUTHORS, AFFILIATIONS, ORCIDs, CORRESPONDING AUTHOR — added before submission]]')
    body = md[md.index('---') + 1:]
    # citations renumbered by first appearance (body first, then the table files); the .md sources keep the verified-list numbers
    refs = references(); rn = Renumberer(refs); body = rn.convert('\n'.join(body)).split('\n')
    table_names = ('Table1_software.md', 'Table2_sequencing_assembly.md', 'Table3_annotation.md') if src.endswith('NOTE.md') else ('Table2_sequencing_assembly.md', 'Table3_annotation.md')
    table_text = {}
    for name in table_names:
        f = os.path.join(MAN, 'tables', name)
        if os.path.exists(f):
            t = open(f, encoding='utf-8').read()
            for a, b in (table_subs or {}).items(): t = t.replace(a, b)
            table_text[name] = rn.convert(t)
    i = 0; in_refs = False; in_legends = False; used = set()
    while i < len(body):
        l = body[i]
        if l.startswith('## References'):
            in_refs = True; i += 1
            while i < len(body) and not body[i].startswith('## '): i += 1
            continue
        if l.startswith('## ') and l[3:].strip() in ('Figure legends', 'Tables'):   # rendered from the files below, not from the notes
            i += 1
            while i < len(body) and not body[i].startswith('## '): i += 1
            continue
        if l.startswith('---'): i += 1; continue
        if l.startswith('## '): doc.add_heading(l[3:].strip(), level=1); i += 1; continue
        if l.startswith('### '): doc.add_heading(l[4:].strip(), level=2); i += 1; continue
        if l.startswith('|'):
            tl = []
            while i < len(body) and body[i].startswith('|'): tl.append(body[i]); i += 1
            md_table(doc, tl); continue
        if l.strip():
            text = l.strip()
            if text.startswith('- '): text = '• ' + text[2:]
            for n in re.findall(r'\[(\d[\d,]*)\]', text): used.update(int(x) for x in n.split(','))
            p = doc.add_paragraph(); rich(p, text)
        i += 1
    # ---- References (only those cited, renumbered in order of first appearance) ----
    doc.add_heading('References', level=1); order = rn.reference_list()
    for n, text in order:
        p = doc.add_paragraph(); p.paragraph_format.line_spacing = 1.15; p.paragraph_format.space_after = Pt(4); rich(p, f'{n}. {text}', 11)
    # ---- Tables ----
    if not extra_tables: table_names = ()
    if table_names: doc.add_page_break(); doc.add_heading('Tables' if src.endswith('NOTE.md') else 'Additional tables (for reference; Table 1 is in the text)', level=1)
    for name in table_names:
        if name not in table_text: continue
        lines = table_text[name].split('\n'); j = 0
        while j < len(lines):
            if lines[j].startswith('**'): p = doc.add_paragraph(); rich(p, lines[j]); j += 1; continue
            if lines[j].startswith('|'):
                tl = []
                while j < len(lines) and lines[j].startswith('|'): tl.append(lines[j]); j += 1
                md_table(doc, tl, 8); continue
            j += 1
    if src.endswith('NOTE.md') and not os.path.exists(os.path.join(MAN, 'tables', 'Table3_annotation.md')): p = doc.add_paragraph(); rich(p, '**Table 3.** [[PENDING: final gene set]]')
    # ---- Figures ----
    doc.add_page_break(); doc.add_heading('Figures', level=1)
    figs = figs if figs is not None else [('Figure 1', 'Fig1_assembly_validation.png', 'Assembly validation (panels complementary to Fig. 1 of the Data Note). (A) GenomeScope 2.0 21-mer spectrum of the HiFi reads (haploid length 2.26 Gb, heterozygosity 0.33 %). (B) Nx curves and cumulative length for Mimp_v3, hap1 and hap2. (C) Merqury QV per chromosome (dotted line: genome-wide QV 67.9). (D) Whole-genome alignment of Mimp_v3 to *Gopherus evgoodei* rGopEvg1 (minimap2 asm20, alignments ≥ 50 kb, MAPQ ≥ 30).'),
            ('Figure 2', 'Fig2_annotation_haplotypes_mitogenome.png', 'Annotation, haplotypes and mitogenome. (A) Protein-coding genes per chromosome (bars) and gene density (line). (B) SNV density per chromosome between each haplotype assembly and the primary assembly. (C) Mitochondrial genome (18,133 bp, circular; 13 protein-coding genes, 22 tRNAs, 2 rRNAs); the tRNA-Thr/tRNA-Pro tandem array lies between CYTB and the control region. (D) Length heteroplasmy of that array: net insertion/deletion length within the array for each of 1,529 spanning HiFi reads (0 bp = four-unit consensus, +300 bp = five units).'),
            ('Supplementary Figure S1', 'FigS1_haplotype_hic_maps.png', 'Hi-C contact maps of the two haplotype assemblies (PretextMap).')]
    for name, f, cap in figs:
        path = f if (f and os.path.isabs(f)) else (os.path.join(FIG, f) if f else None)
        if path and os.path.exists(path): doc.add_picture(path, width=Inches(6.0))
        p = doc.add_paragraph(style='Caption2'); rich(p, f'**{name}.** ' + cap, 10)
    doc.save(out); print('wrote', out, '| cited references:', len(order))

BMC_FIGS = [('Fig. 1', os.path.join(BMCFIG, 'Fig1_assembly_overview.png'), 'The Mimp_v3 assembly. (a) Hi-C contact map of the release (chr1–chr26, then unplaced scaffolds). (b) The 26 chromosomes with contig junctions, telomeric arrays and their *Gopherus evgoodei* homolog. (c) Repeat landscape: genome fraction per 1 % divergence bin, by class. (d) Merqury copy-number spectrum (QV 67.9, k-mer completeness 95.1 %).')]
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'bmc':
        main('MANUSCRIPT_DATA_NOTE_BMC.md', os.path.join(MAN, 'Manouria_impressa_genome_v3_BMC_Genomic_Data_DataNote.docx'), BMC_FIGS, extra_tables=False)   # Data Note: Table 1 only (BMC)
    else: main()
