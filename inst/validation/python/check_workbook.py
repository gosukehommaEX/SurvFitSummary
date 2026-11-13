"""Independent check of a workbook written by write_survival_workbook() for
inst/extdata/fit_example.csv (all models of fit_survival_models() with the
defaults), against the data and the reference values of
inst/validation/python/expected (numpy and scipy only).

Usage: python3 -I check_workbook.py <workbook.xlsx> <expected folder> <fit_example.csv>
Prints PASS or FAIL for each check and ALL PASS at the end."""
import csv, math
import numpy as np
import openpyxl
from scipy import stats, integrate

import sys
WORKBOOK, EXPECTED, DATA = sys.argv[1:4]
wb = openpyxl.load_workbook(WORKBOOK, data_only=True)
def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:]]
def expected(name):
    with open(f'{EXPECTED}/{name}.csv') as f:
        out = []
        for r in csv.DictReader(f):
            out.append({k: (float(v) if v not in ('NA', '') and k not in
                            ('fit_type','arm','model','scale','row','col','parameter') else (None if v == 'NA' else v))
                        for k, v in r.items()})
        return out
results = []
def report(label, value, limit):
    ok = value <= limit
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {label}: {value:.3g} (limit {limit:.3g})")

# ---- data
data = {}
with open(DATA) as f:
    for r in csv.DictReader(f):
        data.setdefault(r['TRT01P'], []).append((float(r['AVAL']), 1 - int(r['CNSR'])))
arms = ['Control', 'Treatment']
z = stats.norm.ppf(0.975)

def km(arm):
    t = np.array([a for a, _ in data[arm]]); d = np.array([b for _, b in data[arm]])
    order = np.argsort(t); t, d = t[order], d[order]
    times = np.unique(t[d == 1])
    s, g, out = 1.0, 0.0, []
    for u in times:
        n = np.sum(t >= u); e = np.sum((t == u) & (d == 1))
        s *= 1 - e / n
        if n > e: g += e / (n * (n - e))
        out.append((u, s, g))
    return t, d, out
def km_at(arm, x):
    t, d, out = km(arm)
    s, g = 1.0, 0.0
    for u, su, gu in out:
        if u <= x: s, g = su, gu
    return s, s * math.sqrt(g)

# Data summary
ds = {r['ARM']: r for r in sheet('Data summary')}
diff = 0
for a in arms:
    t = [x for x, _ in data[a]]; d = [y for _, y in data[a]]
    diff = max(diff, abs(ds[a]['n'] - len(t)), abs(ds[a]['events'] - sum(d)),
               abs(ds[a]['censored'] - (len(t) - sum(d))))
    _, _, out = km(a)
    med = min(u for u, s, _ in out if s <= 0.5)
    diff = max(diff, abs(ds[a]['median'] - med))
report('Data summary: counts and KM median (abs diff)', diff, 1e-9)

# KM table
diff = 0; n_rows = 0
for r in sheet('KM table'):
    a, x = r['ARM'], r['time']
    t = np.array([u for u, _ in data[a]])
    s, se = km_at(a, x)
    lo, hi = s * math.exp(-z * se / s), min(1, s * math.exp(z * se / s))
    diff = max(diff, abs(r['n_risk'] - np.sum(t >= x)), abs(r['surv'] - s),
               abs(r['se'] - se), abs(r['lcl'] - lo), abs(r['ucl'] - hi))
    n_rows += 1
report(f'KM table ({n_rows} rows): n at risk, S, Greenwood SE, log CI', diff, 1e-9)

# Cox model (no ties): Newton-Raphson on the partial likelihood
t = np.array([u for a in arms for u, _ in data[a]])
d = np.array([e for a in arms for _, e in data[a]])
xx = np.array([0.0] * len(data['Control']) + [1.0] * len(data['Treatment']))
b = 0.0
for _ in range(50):
    u_sum, info = 0.0, 0.0
    for i in np.where(d == 1)[0]:
        risk = t >= t[i]; w = np.exp(b * xx[risk])
        m = np.sum(w * xx[risk]) / np.sum(w)
        u_sum += xx[i] - m
        info += np.sum(w * xx[risk] ** 2) / np.sum(w) - m ** 2
    b += u_sum / info
se_b = 1 / math.sqrt(info)
hr = sheet('Hazard ratios')[0]
diff = max(abs(hr['hr'] - math.exp(b)) / math.exp(b),
           abs(hr['lcl'] - math.exp(b - z * se_b)) / math.exp(b),
           abs(hr['ucl'] - math.exp(b + z * se_b)) / math.exp(b))
report('Hazard ratio and Wald CI (relative diff)', diff, 1e-6)

# Fit statistics
gof = {(r['fit_type'], r['arm'], r['model']): r for r in sheet('Fit statistics')}
diff = 0; n = 0
for e in expected('gof'):
    k = (e['fit_type'], e['arm'], e['model'])
    if k in gof:
        n += 1
        diff = max(diff, abs(gof[k]['loglik'] - e['loglik']), abs(gof[k]['AIC'] - e['AIC']),
                   abs(gof[k]['BIC'] - e['BIC']))
report(f'Fit statistics ({n} models in the reference): loglik, AIC, BIC (abs diff)', diff, 1e-3)

# Parameters
for scale, name in [('natural', 'Parameters'), ('transformed', 'Parameters (transformed)')]:
    par = {(r['fit_type'], r['arm'], r['model'], r['parameter']): r for r in sheet(name)}
    de, dse, n = 0, 0, 0
    for e in expected('parameters'):
        if e['scale'] != scale: continue
        k = (e['fit_type'], e['arm'], e['model'], e['parameter'])
        if k in par:
            n += 1
            de = max(de, abs(par[k]['est'] - e['est']) / e['se'])
            dse = max(dse, abs(par[k]['se'] - e['se']) / e['se'])
    report(f'{name} ({n} parameters): |est - ref| / se', de, 0.02)
    report(f'{name}: relative diff of se', dse, 5e-3 * 10)
    # covariance columns against the reference covariance
    rows = sheet(name)
    groups = {}
    for i, r in enumerate(rows):
        groups.setdefault((r['fit_type'], r['arm'], r['model']), []).append(r)
    ref = {(e['fit_type'], e['arm'], e['model'], e['row'], e['col']): e['value']
           for e in expected('covariance') if e['scale'] == scale}
    dc, ddiag, n = 0, 0, 0
    for k, rs in groups.items():
        for i, ri in enumerate(rs):
            ddiag = max(ddiag, abs(ri[f'cov_{i + 1}'] - ri['se'] ** 2) / ri['se'] ** 2)
            for j, rj in enumerate(rs):
                key = k + (ri['parameter'], rj['parameter'])
                if key in ref:
                    n += 1
                    scale_ij = math.sqrt(abs(ref[k + (ri['parameter'], ri['parameter'])] *
                                             ref[k + (rj['parameter'], rj['parameter'])]))
                    dc = max(dc, abs(ri[f'cov_{j + 1}'] - ref[key]) / scale_ij)
            for j in range(len(rs), 6):
                if ri[f'cov_{j + 1}'] is not None: dc = float('inf')
    report(f'{name}: covariance columns vs reference ({n} entries, / sd_i sd_j)', dc, 0.06)
    report(f'{name}: variance on the diagonal = se^2 (relative)', ddiag, 1e-10)

# Covariance and Cholesky sheets: L lower triangular and L L' = covariance
cov, chol = {}, {}
for r in sheet('Covariance'):
    cov.setdefault((r['fit_type'], r['arm'], r['model'], r['scale']), {})[(r['row'], r['col'])] = r['value']
for r in sheet('Cholesky'):
    chol.setdefault((r['fit_type'], r['arm'], r['model'], r['scale']), {})[(r['row'], r['col'])] = r['value']
dl, du, n = 0, 0, 0
for k, c in cov.items():
    names = list(dict.fromkeys(rc[0] for rc in c))
    S = np.array([[c[(a, b2)] for b2 in names] for a in names], dtype=float)
    if np.isnan(S).any(): continue
    L = np.array([[chol[k][(a, b2)] for b2 in names] for a in names], dtype=float)
    n += 1
    du = max(du, np.max(np.abs(np.triu(L, 1))))
    dl = max(dl, np.max(np.abs(L @ L.T - S)) / np.max(np.abs(np.diag(S))))
report(f'Cholesky ({n} matrices): upper triangle', du, 0)
report('Cholesky: L L\' - covariance (relative to the largest variance)', dl, 1e-10)

# Landmark survival against the reference predictions
pred = {(r['fit_type'], r['arm'], r['model'], r['time']): r for r in sheet('Landmark survival')}
ds_, dk, n = 0, 0, 0
for e in expected('predictions'):
    k = (e['fit_type'], e['arm'], e['model'], e['time'])
    if k in pred:
        n += 1
        ds_ = max(ds_, abs(pred[k]['est'] - e['surv']))
        if e['km'] is not None: dk = max(dk, abs(pred[k]['km'] - e['km']))
        elif pred[k]['km'] is not None: dk = float('inf')
report(f'Landmark survival ({n} predictions): S(t) vs reference (abs)', ds_, 1e-3)
report('Landmark survival: KM column vs reference (abs)', dk, 1e-9)

# Survival summary: median and mean against the reference; RMST to 120
# months of the exponential and Weibull fits by numerical integration
ss = {(r['fit_type'], r['arm'], r['model']): r for r in sheet('Survival summary')}
dm, n = 0, 0
for e in expected('summary'):
    k = (e['fit_type'], e['arm'], e['model'])
    if k in ss:
        n += 1
        dm = max(dm, abs(ss[k]['median'] - e['median']) / e['median'])
        if e['mean'] is not None and ss[k]['mean'] is not None:
            dm = max(dm, abs(ss[k]['mean'] - e['mean']) / e['mean'])
report(f'Survival summary ({n} models): median and mean (relative)', dm, 2e-3)
par = {(r['fit_type'], r['arm'], r['model'], r['parameter']): r['est'] for r in sheet('Parameters')}
dr = 0
for a in arms:
    lam = par[('independent', a, 'Exponential', 'rate')]
    rm = integrate.quad(lambda u: math.exp(-lam * u), 0, 120)[0]
    dr = max(dr, abs(ss[('independent', a, 'Exponential')]['rmst'] - rm) / rm)
    sh, sc = par[('independent', a, 'Weibull', 'shape')], par[('independent', a, 'Weibull', 'scale')]
    rm = integrate.quad(lambda u: math.exp(-(u / sc) ** sh), 0, 120)[0]
    dr = max(dr, abs(ss[('independent', a, 'Weibull')]['rmst'] - rm) / rm)
report('RMST to 120 months, exponential and Weibull (relative)', dr, 1e-6)

# Knots
kn = {(r['fit_type'], r['arm'], r['n_knots'], r['knot']): r for r in sheet('Spline knots')}
dk, n = 0, 0
for e in expected('knots'):
    k = (e['fit_type'], e['arm'], e['n_knots'], e['knot'])
    if k in kn:
        n += 1
        dk = max(dk, abs(kn[k]['log_time'] - e['log_time']), abs(kn[k]['n_risk'] - e['n_risk']),
                 abs(kn[k]['surv_est'] - e['surv_est']), abs(kn[k]['surv_se'] - e['surv_se']))
report(f'Spline knots ({n} knots): log time, n at risk, KM, SE (abs)', dk, 1e-6)

# Print-ready names and model status
names = sorted(set(r['model'] for r in sheet('Fit statistics')))
print('Model names:', names)
bad = [r for r in sheet('Model status') if r['status'] != 'ok']
print('Models not ok:', [(r['fit_type'], r['arm'], r['model'], r['status'], r['message']) for r in bad])
print('ALL PASS' if all(results) else 'SOME FAIL')
