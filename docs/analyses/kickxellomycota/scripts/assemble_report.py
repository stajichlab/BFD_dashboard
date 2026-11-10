#!/usr/bin/env python3
"""Assemble the report from the two text parts and the table fragments."""
import re
t=open('kx_scripts/report_part1.txt').read()+open('kx_scripts/report_part2.txt').read()
t=re.sub(r'\{\{(frag_\w+)\}\}',lambda m:open(f'kx_data/{m.group(1)}.md').read().rstrip(),t)
open('fable_A_kickxello_and_agenda.md','w').write(t); print(len(t),'chars;', t.count('![') ,'figures;', len(re.findall(r'\{\{',t)),'unfilled')
