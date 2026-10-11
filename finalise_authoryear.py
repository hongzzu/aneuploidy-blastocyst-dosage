# -*- coding: utf-8 -*-
"""Finalise the citation conversion: handle the 3 ambiguous cases, sort references,
renumber in-text as author-year, write the v6 manuscript."""
import io, re

P = r'F:\xinlunwen\05_Manuscript\论文2_终稿_v3.md'
t = io.open(P, encoding='utf-8').read()

L = {1:('Hassold',2001),2:('Nagaoka',2012),3:('Franasiak',2014),4:('Petropoulos',2016),
     5:('Blakeley',2015),6:('Meistermann',2021),7:('Wang',2024),8:('Squair',2021),
     9:('Zimmerman',2021),10:('Murphy',2022),11:('Carioscia',2026),12:('Sun',2023),
     13:('GTEx Consortium',2020),14:('Võsa',2021),15:('Liu',2019),16:('Benjamini',1995),
     17:('Zöllner',2007),18:('Burgess',2013),19:('Burgess',2011),
     20:('Giambartolomei',2014),21:('Subramanian',2005),22:('Fang',2023),
     23:('Torres',2007),24:('Williams',2008),25:('Santaguida',2015),26:('Stingele',2012),
     27:('Skrivankova',2021),28:('Skrivankova',2021),29:('Starostik',2020),
     30:('Wallace',2021),31:('Iwamoto-Stohl',2026),32:('Steinthorsdottir',2024),
     33:('Brion',2013),34:('Gupta',2024),35:('Alfeghaly',2025),36:('Wei',2025)}
L[27] = ('Skrivankova', '2021a')
L[28] = ('Skrivankova', '2021b')

def lab(nums):
    out = []
    for n in nums:
        s, y = L[n]
        out.append('%s, %s' % (s, y) if n == 13 else '%s et al., %s' % (s, y))
    seen, res = set(), []
    for x in out:
        if x not in seen:
            seen.add(x); res.append(x)
    return '; '.join(res)

def expand(spec):
    out = []
    for part in spec.split(','):
        part = part.strip()
        if '–' in part:
            a, b = part.split('–'); out.extend(range(int(a), int(b)+1))
        elif '-' in part:
            a, b = part.split('-'); out.extend(range(int(a), int(b)+1))
        else:
            out.append(int(part))
    return out

body, refs = t.split('## References')

# 1. ambiguous (27,28) STROBE-MR -> author-year
o = u'**STROBE-MR** (27,28)'
n = u'**STROBE-MR** (%s)' % lab([27, 28])
assert t.count(o) == 1, t.count(o)
body = body.replace(o, n)

# 2. keep (1, 7, 13, 16, 22) and (3–6) as literal numbers (chromosome counts, not citations)

# 3. known citation specs
CIT = ['1,2','3','4–6,29','7','8–10','11','11,12','12','15','16','17','18','19','20,30',
       '21,22','23–26','31','32','33','34','35','36']
for spec in CIT:
    nums = expand(spec)
    old = '(%s)' % spec
    if body.count(old) == 0:
        continue
    body = body.replace(old, '(%s)' % lab(nums))

io.open(P, 'w', encoding='utf-8').write(body + '## References' + refs)
print('in-text citations converted to author-year')

# 4. re-sort the reference list alphabetically and re-emit
t = io.open(P, encoding='utf-8').read()
body, refs = t.split('## References')
entries = {}
for m in re.finditer(r'^(\d+)\.\s+(.*?)(?=\n\d+\.\s|\n## |\Z)', refs, flags=re.S | re.M):
    entries[int(m.group(1))] = ' '.join(m.group(2).split())

def sortkey(num):
    txt = entries[num]
    m = re.match(r"([A-Za-zÀ-ÿ'’\-]+)", txt)
    sur = m.group(1) if m else ''
    sur = sur.lower().replace('ã', 'a').replace('ö', 'o').replace('õ', 'o')
    if num == 13: sur = 'gtex consortium'
    ym = re.search(r'\b(19|20)\d{2}\b', txt)
    yr = ym.group(0) if ym else '0000'
    if num == 27: yr = '2021a'
    if num == 28: yr = '2021b'
    return (sur, yr, num)

order = sorted(entries.keys(), key=sortkey)
lines = []
prev = None
for num in order:
    txt = entries[num]
    m = re.match(r"([A-Za-zÀ-ÿ'’\-]+)", txt)
    sur = m.group(1) if m else ''
    yr = re.search(r'\b(19|20)\d{2}\b', txt)
    yr = yr.group(0) if yr else ''
    lines.append(u'%s %s. %s' % (sur, yr, txt))
new_refs = u'## References\n\n' + u'\n\n'.join(lines) + u'\n'
t = body + new_refs
io.open(P, 'w', encoding='utf-8').write(t)
print('references re-sorted alphabetically: %d entries' % len(order))
