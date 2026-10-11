# -*- coding: utf-8 -*-
"""Merge Table1_data_sources.docx + Table2_top_dosage_genes.docx into Tables.docx."""
import io
import re
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

M = r'F:\xinlunwen\05_Manuscript'
OUT = M + '\\Tables.docx'

# read the source markdown tables (single source of truth) rather than merging docx XML
t1 = io.open(M + '\\Table1_data_sources.md', encoding='utf-8').read()
t2 = io.open(M + '\\Table2_top_dosage_genes.md', encoding='utf-8').read()

doc = Document()
style = doc.styles['Normal']
style.font.name = 'Times New Roman'
style.font.size = Pt(10)

for sec in doc.sections:
    # landscape gives the wide tables room
    from docx.enum.section import WD_ORIENT
    sec.orient = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    sec.left_margin = Inches(0.6)
    sec.right_margin = Inches(0.6)
    sec.top_margin = Inches(0.7)
    sec.bottom_margin = Inches(0.7)

BOLD = re.compile(r'\*\*(.+?)\*\*')
ITAL = re.compile(r'(?<!\*)\*([^*\n]+?)\*(?!\*)')
CODE = re.compile(r'`([^`]+)`')


def add_runs(par, text, base_bold=False, size=None):
    tokens = re.split(r'(\*\*.+?\*\*|`[^`]+`|(?<!\*)\*[^*\n]+?\*(?!\*))', text)
    for tok in tokens:
        if not tok:
            continue
        m = BOLD.fullmatch(tok)
        if m:
            r = par.add_run(m.group(1)); r.bold = True
        elif CODE.fullmatch(tok):
            r = par.add_run(m.group(1) if (m := CODE.fullmatch(tok)) else tok)
            r.font.name = 'Courier New'
        elif ITAL.fullmatch(tok):
            r = par.add_run(ITAL.fullmatch(tok).group(1)); r.italic = True
        else:
            r = par.add_run(tok)
            if base_bold:
                r.bold = True
        if size:
            r.font.size = Pt(size)


def emit_table(lines):
    """Render a markdown table block as a Word table."""
    rows = []
    for l in lines:
        s = l.strip()
        if s.startswith('|') and not re.match(r'^\|[\s:\-|]+\|$', s):
            cells = [c.strip() for c in s.strip('|').split('|')]
            rows.append(cells)
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    tb = doc.add_table(rows=len(rows), cols=ncol)
    tb.style = 'Table Grid'
    for ri, row in enumerate(rows):
        for ci in range(ncol):
            cell = tb.cell(ri, ci)
            txt = row[ci] if ci < len(row) else ''
            par = cell.paragraphs[0]
            par.paragraph_format.space_after = Pt(2)
            add_runs(par, txt, base_bold=(ri == 0), size=8)
    doc.add_paragraph()


def emit_block(md, title):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    r = p.add_run(title)
    r.bold = True; r.font.size = Pt(13)
    p.paragraph_format.space_after = Pt(10)

    lines = md.split('\n')
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1; continue
        if s.startswith('|'):
            block = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                block.append(lines[i]); i += 1
            emit_table(block)
        elif s.startswith('#'):
            p = doc.add_paragraph()
            r = p.add_run(s.lstrip('# ').strip())
            r.bold = True; r.font.size = Pt(11)
        elif s.startswith('**'):
            p = doc.add_paragraph()
            add_runs(p, s)
        else:
            p = doc.add_paragraph()
            add_runs(p, s)
        i += 1


emit_block(t1, 'Table 1. Data sources, accessions and analytical role')
doc.add_page_break()
emit_block(t2, 'Table 2. The 50 dosage-sensitive genes with the strongest evidence in the cleaned set')

doc.save(OUT)
print('saved', OUT)
d2 = Document(OUT)
print('tables in merged doc:', len(d2.tables))
print('paragraphs:', len(d2.paragraphs))
