# -*- coding: utf-8 -*-
"""Check the rendered DOCX files for leftover Markdown markers."""
import re
import sys
from docx import Document

M = r'F:\xinlunwen\05_Manuscript'
FILES = ['\u8bba\u65872_\u7ec8\u7a3f_v3', 'Supplementary_Materials',
         'Table1_data_sources', 'Table2_top_dosage_genes']

# a marker is: a run of **, a lone * adjacent to a word boundary, a backtick, or ~~
MARK = re.compile(r'\*\*|\*|`|~~')


def texts_of(path):
    d = Document(path)
    out = [p.text for p in d.paragraphs]
    for tb in d.tables:
        for row in tb.rows:
            for c in row.cells:
                out.append(c.text)
    return out


for f in FILES:
    ts = texts_of(M + '\\' + f + '.docx')
    bad = [s for s in ts if MARK.search(s)]
    print('%-28s cells/paragraphs containing a marker: %d' % (f, len(bad)))
    for b in bad[:8]:
        print('      %s' % b[:120])
    print()
