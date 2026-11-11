import os  # AO_BASE: directory holding ao/, fable_C_ao_results/, figures
"""Assemble the markdown report from prose fragments (report/*.txt) and rendered tables (report/tables.json, report/jk.md)."""
import json, re
T = json.load(open('report/tables.json')); T['jk'] = open('report/jk.md').read()
main = open('report/01_main.txt').read(); meth = open('report/02_methods.txt').read()
i = main.index('## 3. Genome size'); doc = main[:i] + meth + '\n' + main[i:]
missing = [k for k in re.findall(r'\{\{(\w+)\}\}', doc) if k not in T]; assert not missing, missing
doc = re.sub(r'\{\{(\w+)\}\}', lambda m: T[m.group(1)], doc)
out = os.environ.get('AO_BASE', '.') + '/fable_C_antarctolithica.md'
open(out, 'w').write(doc); print(len(doc), 'chars ->', out)
