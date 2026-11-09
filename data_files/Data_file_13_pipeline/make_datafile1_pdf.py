#!/usr/bin/env python3
"""Data file 1 PDF: long-form Markdown (manuscript/MANUSCRIPT_DATA_NOTE.md) + Tables 1-3 + long-form figures → HTML → headless Chrome PDF.
Only cited references are listed, renumbered by first appearance (refnumber.py). Requires Google Chrome on macOS; python-markdown."""
import markdown, re, base64, subprocess, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from refnumber import Renumberer
BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); os.chdir(BASE)
MAN = 'docs/genome_paper_new_submission/manuscript'; FIG = 'docs/genome_paper_new_submission/figures/longform'
md = open(f'{MAN}/MANUSCRIPT_DATA_NOTE.md').read().split('---', 1)[1]
refs = {int(m.group(1)): m.group(2) for m in re.finditer(r'^(\d+)\.\s+(.*)', open(f'{MAN}/references_verified.md').read(), re.M)}
rn = Renumberer(refs); body, _, rest = md.partition('## References'); body = rn.convert(body)
tables = '\n\n## Tables\n\n' + rn.convert('\n\n'.join(open(f'{MAN}/tables/{t}.md').read() for t in ['Table1_software', 'Table2_sequencing_assembly', 'Table3_annotation']))
cited = rn.reference_list()
md = body + '## References\n\n' + '\n'.join(f'{n}. {t}' for n, t in cited) + '\n' + re.sub(r'^.*?(?=\n## |\Z)', '', rest, count=1, flags=re.S)
img = lambda p: f'<img src="data:image/png;base64,{base64.b64encode(open(p, "rb").read()).decode()}" style="max-width:100%">'
figs = '\n\n## Figures\n\n' + '\n\n'.join(f'<p><b>{t}</b></p>{img(f"{FIG}/{f}")}' for t, f in [('Figure 1', 'Fig1_assembly_validation.png'), ('Figure 2', 'Fig2_annotation_haplotypes_mitogenome.png'), ('Supplementary Figure S1', 'FigS1_haplotype_hic_maps.png')])
body = markdown.markdown(md + tables + figs, extensions=['tables'])
html = ('<html><head><meta charset="utf-8"><title>Data file 1</title><style>body{font-family:Georgia,serif;font-size:11pt;line-height:1.45;max-width:7in;margin:0.6in auto} '
        'table{border-collapse:collapse;font-size:8.5pt;width:100%} td,th{border:1px solid #999;padding:3px 5px;vertical-align:top} code{font-size:9pt} h1{font-size:16pt} h2{font-size:13pt;margin-top:1.2em} h3{font-size:11.5pt}</style></head><body>'
        '<p style="color:#666;font-size:9pt">Data file 1 of: A chromosome-level, haplotype-resolved genome assembly of the impressed tortoise, <i>Manouria impressa</i> (BMC Genomic Data). Extended methods, validation and supplementary tables. Items in [[ ]] are pending accessions or author-supplied text.</p>' + body + '</body></html>')
tmp = f'{MAN}/.datafile1.html'; open(tmp, 'w').write(html)
out = f'{BASE}/{MAN}/Data_file_01_extended_description.pdf'
subprocess.run(['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '--headless=new', '--disable-gpu', '--no-pdf-header-footer', f'--print-to-pdf={out}', f'file://{BASE}/{tmp}'], capture_output=True)
os.remove(tmp); print('wrote', out, os.path.getsize(out), 'bytes; references listed:', len(cited))
