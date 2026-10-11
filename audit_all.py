# -*- coding: utf-8 -*-
"""Comprehensive consistency audit of the manuscript and its companions."""
import io
import json
import os
import re

M = r'F:\xinlunwen\05_Manuscript'
RES = r'F:\xinlunwen\results'
V3 = io.open(os.path.join(M, u'\u8bba\u65872_\u7ec8\u7a3f_v3.md'), encoding='utf-8').read()
CL = io.open(os.path.join(M, 'Cover_Letter_v2.md'), encoding='utf-8').read()
HL = io.open(os.path.join(M, 'Highlights_v2.md'), encoding='utf-8').read()
CK = io.open(os.path.join(M, 'Submission_Checklist.md'), encoding='utf-8').read()
SP = io.open(os.path.join(M, 'Supplementary_Materials.md'), encoding='utf-8').read()
T1 = io.open(os.path.join(M, 'Table1_data_sources.md'), encoding='utf-8').read()
T2 = io.open(os.path.join(M, 'Table2_top_dosage_genes.md'), encoding='utf-8').read()

print('=' * 78)
print('1. WORD COUNTS')
print('=' * 78)
def prose(t, a, b):
    s = t.split(a)[1].split(b)[0]
    s = '\n'.join(l for l in s.split('\n') if not l.strip().startswith('|'))
    return len(s.split())
ab = len(V3.split('## Abstract')[1].split('**Keywords**')[0].split())
i = prose(V3, '## 1. Introduction', '## 2. Methods')
m = prose(V3, '## 2. Methods', '## 3. Results')
r = prose(V3, '## 3. Results', '## 4. Discussion')
d = prose(V3, '## 4. Discussion', '## 5. Declarations')
ric = 0
print('abstract %d | intro %d | methods %d | results %d | discussion %d | TOTAL %d'
      % (ab, i, m, r, d, i + m + r + d))
for label, tgt, val in [('abstract', 250, ab), ('methods', 1200, m),
                        ('results', 1600, r), ('discussion', 1300, d),
                        ('main total', 5000, i + m + r + d)]:
    print('   %-12s %5d / %-5d %s' % (label, val, tgt, 'OK' if val <= tgt else 'OVER'))

print()
print('=' * 78)
print('2. FIGURE CITATIONS vs PANELS AVAILABLE')
print('=' * 78)
legs = {}
for n in range(1, 6):
    seg = V3.split('**Figure %d.' % n)[1]
    nxt = re.search(r'\*\*Figure \d\.', seg)
    legs[n] = seg[:nxt.start()] if nxt else seg
avail = {n: sorted(set(re.findall(r'\*\*\(([a-f])\)\*\*', legs[n]))) for n in legs}
body = V3.split('## Figure legends')[0]
cited = {}
for n in range(1, 6):
    pat = re.compile(r'Fig(?:ure)?\.?\s*%d\s*([a-f]?(?:\s*,\s*[a-f])*)' % n)
    letters = set()
    for mm in pat.finditer(body):
        letters |= set(re.findall(r'[a-f]', mm.group(1)))
    cited[n] = sorted(letters)
for n in range(1, 6):
    have = set(avail[n])
    used = set(cited[n])
    orphan = sorted(have - used)
    ghost = sorted(used - have)
    print('Figure %d: panels %s | cited %s %s%s'
          % (n, ''.join(avail[n]), ''.join(cited[n]) or '(whole figure)',
             ('  ORPHAN PANELS: ' + ','.join(orphan)) if orphan else '',
             ('  !! CITES NON-EXISTENT: ' + ','.join(ghost)) if ghost else ''))

print()
print('=' * 78)
print('3. TABLE CITATIONS')
print('=' * 78)
for t in ['Table 1', 'Table 2', 'Table S1', 'Table S2', 'Table S3', 'Table S4',
          'Table S5', 'Table S6']:
    print('  %-10s cited in body: %d' % (t, len(re.findall(re.escape(t) + r'\b', body))))
print('  Table 2 rows in source:', T2.count('\n| ') - 1)

print()
print('=' * 78)
print('4. KEY NUMBERS vs RESULT FILES')
print('=' * 78)
ms = json.load(io.open(os.path.join(RES, 'manuscript_stats.json'), encoding='utf-8'))
cs = json.load(io.open(os.path.join(RES, 'concordance_stats.json'), encoding='utf-8'))
checks = [
    ('198', '198'), ('14,082', '14,082'), ('2,833', '2,833'), ('3,635', '3,635'),
    ('40,267', '40,267'), ('6,079', '6,079'), ('15.1%', '15.1%'),
    ('+0.479', '+0.479'), ('+0.186', '+0.186'), ('+0.772', '+0.772'),
    ('0.0082', '0.0082'), ('4,310', '4,310'), ('98.4%', '98.4%'),
    ('94.9%', '94.9%'), ('+0.620', '+0.620'), ('1.09%', '1.09%'),
    ('103,906', '103,906'), ('9,541,215', '9,541,215'),
    ('15,333,206', '15,333,206'), ('16,685', '16,685'), ('19,030', '19,030'),
    ('479', '479'), ('1,470', '1,470'),
]
for a, b in checks:
    print('  %-12s in body: %-3d  in abstract: %d'
          % (a, V3.count(a), V3.split('## Abstract')[1].split('**Keywords**')[0].count(a)))

print()
print('=' * 78)
print('5. CROSS-DOCUMENT CLAIMS')
print('=' * 78)
print('cover letter "main text 4,900 words"      :', '4,900' in CL)
print('cover letter "abstract 233 words"         :', '233' in CL)
print('cover letter "5 main figures"             :', '5 main figures' in CL)
print('cover letter ORCID                        :', '0009-0005-3900-9648' in CL)
print('code statement in manuscript              :',
      'available from the corresponding author on reasonable request' in V3)
print('code statement in cover letter            :',
      'available from the corresponding author on reasonable request' in CL)
print('code statement in Table 1                 :',
      'available from the corresponding author on reasonable request' in T1)
print('code statement in Supplementary           :',
      'available from the corresponding author on reasonable request' in SP)
print('supplementary cites Figure S1/S2/S3       :',
      len(re.findall(r'Fig(?:ure)?\s*S\d', SP)))
print('supplementary cites a Figure at all       :',
      len(re.findall(r'Fig(?:ure)?\s*\d', SP)))

print()
print('=' * 78)
print('6. MARKDOWN / STRUCTURE SANITY')
print('=' * 78)
print('unclosed bold (**) count (should be even) :', V3.count('**') % 2 == 0, V3.count('**'))
print('duplicate consecutive sentences           :', end=' ')
ss = [x for x in re.split(r'(?<=[.])\s+', body) if len(x) > 60]
dups = [s for s in set(ss) if ss.count(s) > 1]
print(len(dups))
for x in dups[:5]:
    print('     DUP: %s' % x[:100])
print('placeholders [To be completed.]           :', V3.count('[To be completed.]'))
print('CJK characters in manuscript              :',
      len([c for c in V3 if u'\u4e00' <= c <= u'\u9fff']))
print('references numbered 1..30                 :',
      [int(x) for x in re.findall(r'^(\d+)\. ', V3, flags=re.M)][-30:] == list(range(1, 31)))
