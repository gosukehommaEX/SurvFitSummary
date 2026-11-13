"""Independent check of a report written by write_survival_report() for
inst/extdata/fit_example.csv (all models of fit_survival_models() with the
defaults): numbers in the tables against the data (own Kaplan-Meier and
Cox) and the reference values of inst/validation/python/expected,
numbering of the tables and figures.

Usage: python3 -I check_report.py <report.docx> <expected folder> <fit_example.csv>
         [--csv <verdicts.csv>] [--mutate]
Prints PASS or FAIL for each check and ALL PASS at the end; --csv also writes
the verdicts (with the columns of verdict_table() in verdict.R). With
--mutate, every number in the cells of the tables, except in the columns that
identify a row, is shifted by 0.01 * max(1, |value|) and printed with its
number of decimals (at least one unit of its last digit), as in the mutation
run of run_validation.R, and the checks of the numbering of the tables and
figures, which involve no number of the tables, are left out: every check
must then FAIL, which shows that the checks are not vacuous. A check that
compares no row FAILs."""
import argparse, csv, math, re
import numpy as np
import docx
from docx.oxml.ns import qn
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('report'); parser.add_argument('expected'); parser.add_argument('data')
parser.add_argument('--csv'); parser.add_argument('--mutate', action='store_true')
args = parser.parse_args()
path, EXPECTED, DATA = args.report, args.expected, args.data
# Columns that identify a row: not shifted by --mutate
KEYS = {'Treatment', 'Time (months)', 'Distribution', 'Model', 'n_knots', 'knot'}
def shift_number(m):
    text = m.group(0)
    decimals = len(text.split('.')[1]) if '.' in text else 0
    v = float(text)
    new = f'{v + 0.01 * max(1.0, abs(v)):.{decimals}f}'
    if new == text:
        new = f'{v + 10 ** -decimals:.{decimals}f}'
    return new
def shifted(column, text):
    if not args.mutate or column in KEYS:
        return text
    return re.sub(r'-?\d+(?:\.\d+)?', shift_number, text)
d = docx.Document(path)
paras, tables, images, caption = [], [], 0, None
ti = 0
for child in d.element.body.iterchildren():
    tag = child.tag.split('}')[1]
    if tag == 'p':
        text = docx.text.paragraph.Paragraph(child, d).text
        images += len(child.findall('.//' + qn('w:drawing')))
        if text: paras.append(text)
        if text.startswith('Table '): caption = text
    elif tag == 'tbl':
        t = d.tables[ti]; ti += 1
        cells = [[c.text for c in r.cells] for r in t.rows]
        cells = cells[:1] + [[shifted(h, c) for h, c in zip(cells[0], r)] for r in cells[1:]]
        tables.append((caption, cells))
results = []
def report(label, ok, detail='', numbering=False):
    # numbering: a check of the numbering, left out with --mutate
    if numbering and args.mutate:
        return
    results.append((label, detail, ok)); print(f"{'PASS' if ok else 'FAIL'}  {label} {detail}")
def table(name):
    for cap, rows in tables:
        if cap is not None and re.sub(r'^Table \d+: ', '', cap) == name:
            return [dict(zip(rows[0], r)) for r in rows[1:]]
    raise KeyError(name)
def expected(name):
    with open(f'{EXPECTED}/{name}.csv') as f:
        return list(csv.DictReader(f))

# numbering
tab_no = [int(m.group(1)) for p in paras for m in [re.match(r'Table (\d+): ', p)] if m]
fig_no = [int(m.group(1)) for p in paras for m in [re.match(r'Figure (\d+): ', p)] if m]
report('tables numbered 1..n, one caption per table', tab_no == list(range(1, len(tables) + 1)), f'({len(tables)} tables)', numbering=True)
report('figures numbered 1..n, one caption per image', fig_no == list(range(1, images + 1)), f'({images} images)', numbering=True)
# 2 arms, 7 standard models and 9 spline models
expected_figures = 2 * 3 * 2 + 3 + 2 * (2 + 4 * 7) + (2 * 2 + 4 * 7) + 2 * 2 + 2 * 2
report('number of figures', images == expected_figures, f'({images} vs {expected_figures})', numbering=True)
report('no figure failed', not any('could not be drawn' in p for p in paras), numbering=True)

# data and Kaplan-Meier
data = {}
with open(DATA) as f:
    for r in csv.DictReader(f):
        data.setdefault(r['TRT01P'], []).append((float(r['AVAL']), 1 - int(r['CNSR'])))
z = stats.norm.ppf(0.975)
def km_at(arm, x):
    t = np.array([a for a, _ in data[arm]]); e = np.array([b for _, b in data[arm]])
    s, g = 1.0, 0.0
    for u in np.unique(t[(e == 1) & (t <= x)]):
        n = np.sum(t >= u); k = np.sum((t == u) & (e == 1))
        s *= 1 - k / n; g += k / (n * (n - k)) if n > k else 0
    return s, s * math.sqrt(g), np.sum(t >= x)
def km_median(arm):
    t = np.array([a for a, _ in data[arm]]); e = np.array([b for _, b in data[arm]])
    s = 1.0
    for u in np.unique(t[e == 1]):
        s *= 1 - np.sum((t == u) & (e == 1)) / np.sum(t >= u)
        if s <= 0.5: return u
f1 = lambda v: f'{v:.1f}'; f2 = lambda v: f'{v:.2f}'; f3 = lambda v: f'{v:.3f}'
ok = True
for r in table('Data summary'):
    a = r['Treatment']; ev = sum(e for _, e in data[a])
    ok &= r['N patients'] == str(len(data[a])) and r['N events'] == str(ev) \
        and r['N censored'] == str(len(data[a]) - ev)
    ok &= r['Median time to event, months (95% CI)'].split(' ')[0] == f1(km_median(a))
report('Table "Data summary": counts and medians', ok and len(table('Data summary')) == 2)
ok = True; n = 0
for a in ['Control', 'Treatment']:
    for r in table(f'Kaplan-Meier estimates of {a}'):
        x = float(r['Time (months)']); s, se, nr = km_at(a, x); n += 1
        lo, hi = s * math.exp(-z * se / s), min(1, s * math.exp(z * se / s))
        ok &= r['Number at risk'] == str(nr)
        ok &= r['Survival probability (SE)'] == f'{f3(s)} ({f3(se)})'
        ok &= r['Survival probability 95% CI'] == f'{f3(lo)}, {f3(hi)}'
report('Tables "Kaplan-Meier estimates": n at risk, S (SE), CI', ok and n > 0, f'({n} rows)')

# Cox model (no ties)
t = np.array([u for a in data for u, _ in data[a]]); e = np.array([v for a in data for _, v in data[a]])
xx = np.array([1.0 if a == 'Treatment' else 0.0 for a in data for _ in data[a]])
b = 0.0
for _ in range(50):
    U = I = 0.0
    for i in np.where(e == 1)[0]:
        w = np.exp(b * xx[t >= t[i]]); xr = xx[t >= t[i]]; m = np.sum(w * xr) / np.sum(w)
        U += xx[i] - m; I += np.sum(w * xr ** 2) / np.sum(w) - m ** 2
    b += U / I
se = 1 / math.sqrt(I)
r = table('Hazard ratios from Cox proportional hazards models')[0]
p = 2 * stats.norm.sf(abs(b / se))
report('Table "Hazard ratios": HR (CI), p-value',
       r['HR (95% CI)'] == f'{f2(math.exp(b))} ({f2(math.exp(b - z * se))}, {f2(math.exp(b + z * se))})'
       and r['p-value'] == ('<0.001' if p < 0.001 else f3(p)), r['HR (95% CI)'])

# goodness of fit, knots and predictions against the reference
gof = {(g['fit_type'], g['arm'], g['model']): g for g in expected('gof')}
ok = True; n = 0
for arm, fit in [('Control', 'independent'), ('Treatment', 'independent'), ('Both arms', 'dependent')]:
    rows = table(f'Goodness of fit statistics of the standard parametric models ({arm}), in ascending order of AIC')
    aic = [float(r['AIC']) for r in rows]
    ok &= aic == sorted(aic)
    for r in rows:
        g = gof[(fit, arm, r['Distribution'])]; n += 1
        ok &= abs(float(r['AIC']) - float(g['AIC'])) <= 0.05 + 1e-3
        ok &= abs(float(r['BIC']) - float(g['BIC'])) <= 0.05 + 1e-3
report('Tables "Goodness of fit": AIC, BIC and their order', ok and n > 0, f'({n} rows)')
kn = {(k['fit_type'], k['arm'], k['n_knots'], k['knot']): k for k in expected('knots')}
ok = True; n = 0
for arm, fit in [('Control', 'independent'), ('Treatment', 'independent'), ('Both arms', 'dependent')]:
    for r in table(f'Placement of knots ({arm})'):
        k = kn[(fit, arm, r['n_knots'], r['knot'])]; n += 1
        ok &= r['n_risk'] == k['n_risk'] and abs(float(r['time']) - float(k['time'])) <= 0.005 + 1e-6
        ok &= abs(float(r['surv_est']) - float(k['surv_est'])) <= 5e-4 + 1e-9
        ok &= abs(float(r['quantile']) - float(k['quantile'])) <= 5e-4 + 1e-9
report('Tables "Placement of knots"', ok and n > 0, f'({n} rows)')
pred = {}
for p_ in expected('predictions'):
    pred[(p_['fit_type'], p_['arm'], p_['model'], float(p_['time']))] = float(p_['surv'])
summ = {(s_['fit_type'], s_['arm'], s_['model']): s_ for s_ in expected('summary')}
cols = [('Survival at 6 months', 6), ('Survival at 1 year', 12), ('Survival at 2 years', 24),
        ('Survival at 5 years', 60), ('Survival at 10 years', 120)]
ok = True; n = 0
for arm, fit, title in [('Control', 'independent', 'Control'), ('Treatment', 'independent', 'Treatment'),
                        ('Both arms', 'dependent', 'Both arms, Control'), ('Both arms', 'dependent', 'Both arms, Treatment')]:
    for r in table(f'Tests of fit and predicted survival probability ({title})'):
        a = title.split(', ')[-1]; n += 1
        g = gof[(fit, arm, r['Model'])]
        ok &= abs(float(r['AIC']) - float(g['AIC'])) <= 0.05 + 1e-3
        for c, tm in cols:
            ok &= abs(float(r[c]) - pred[(fit, a, r['Model'], tm)]) <= 5e-4 + 1e-3
        s_ = summ[(fit, a, r['Model'])]
        ok &= abs(float(r['Median (months)']) - float(s_['median'])) <= 0.05 + 2e-3 * float(s_['median'])
        ok &= abs(float(r['Mean (months)']) - float(s_['mean'])) <= 0.05 + 2e-3 * float(s_['mean'])
report('Tables "Tests of fit and predicted survival": AIC, median, mean, survival', ok and n > 0, f'({n} rows)')
gt = [p_ for p_ in paras if p_.startswith('The Grambsch-Therneau test')]
print('   G-T text:', gt[0][:90])
oks = [ok for *_, ok in results]
if args.mutate:
    print('MUTATION DETECTED: every check FAILS' if not any(oks) else
          f'MUTATION NOT DETECTED: {sum(oks)} check(s) PASS')
else:
    print('ALL PASS' if all(oks) else 'SOME FAIL')
if args.csv:
    with open(args.csv, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\r\n')
        w.writerow(['script', 'item', 'reference', 'value', 'difference', 'tolerance',
                    'max_difference', 'verdict'])
        for label, detail, ok in results:
            w.writerow(['check_report', f'{label} {detail}'.strip(), 'NA', 'NA', 'NA', 'NA',
                        'NA', 'PASS' if ok else 'FAIL'])
