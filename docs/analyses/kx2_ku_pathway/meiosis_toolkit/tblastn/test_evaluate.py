"""Fixture test for evaluate_tblastn.py: hand-made calls with a known pass and a known fail. Run: python3 -I test_evaluate.py"""
import csv, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ['locustag', 'genus', 'grp', 'gene', 'call']


def run():
    d = tempfile.mkdtemp()
    calls = os.path.join(d, 'calls.tsv')
    rows = []
    # gene G_pass: 4 reference-group positives, 3 recovered (75%) -> PASS ; gene G_fail: 4 positives, 2 recovered (50%) -> FAIL
    for i, c in enumerate(['strong', 'weak', 'strong', 'none']):
        rows.append([f'P{i}', 'g', 'Harpellales', 'G_pass', c])
    for i, c in enumerate(['strong', 'none', 'weak', 'none']):
        rows.append([f'P{i}', 'g', 'Harpellales', 'G_fail', c])
    rows.append(['K1', 'g', 'Kickxellales_core', 'G_pass', 'none'])          # a non-positive genome must not count
    with open(calls, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(COLS)
        w.writerows(rows)
    pos = os.path.join(d, 'pos.tsv')
    with open(pos, 'w') as f:
        f.write('key\tn\tannotated_positive_reference_genomes\nG_pass\t4\tP0,P1,P2,P3\nG_fail\t4\tP0,P1,P2,P3\n')
    subprocess.run([sys.executable, '-I', os.path.join(HERE, 'evaluate_tblastn.py'), calls, pos, os.path.join(d, 'out')], check=True,
                   stdout=subprocess.DEVNULL)
    res = {r['key']: r for r in csv.DictReader(open(os.path.join(d, 'out_control_pass.tsv')), delimiter='\t')}
    assert res['G_pass']['status'].startswith('PASS') and res['G_pass']['pct_recovered'] == '75.0', res['G_pass']
    assert res['G_fail']['status'].startswith('FAIL') and res['G_fail']['pct_recovered'] == '50.0', res['G_fail']


if __name__ == '__main__':
    run()
    print('evaluate_tblastn fixture test passed')
