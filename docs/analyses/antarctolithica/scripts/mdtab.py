"""Tiny helper: DataFrame -> GitHub markdown table."""
import pandas as pd, numpy as np
def md(df, floatfmt=3, index=False):
    d = df.reset_index() if index else df
    def f(v):
        if isinstance(v, (float, np.floating)):
            if np.isnan(v): return ''
            return f'{v:.{floatfmt}g}' if abs(v) < 1e4 else f'{v:,.0f}'
        return str(v)
    out = ['| ' + ' | '.join(map(str, d.columns)) + ' |', '|' + '---|' * len(d.columns)]
    for r in d.itertuples(index=False): out.append('| ' + ' | '.join(f(v) for v in r) + ' |')
    return '\n'.join(out)
