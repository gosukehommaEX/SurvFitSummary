"""Independent reference values for the Kaplan-Meier, Cox and Grambsch-Therneau
tests of SurvFitSummary (tests/testthat/test-summarise_data.R, test-km_table.R,
test-estimate_hazard_ratio.R and test-test_proportional_hazards.R).

The computations use numpy and scipy only (no survival analysis library):
Kaplan-Meier with Greenwood variances, the three pointwise confidence
intervals of survival::survfit(), the median and its confidence limits as the
first time the curve is at or below 0.5, the Cox partial likelihood with the
Efron approximation for ties (Newton-Raphson, optional strata), and the score
test of the Grambsch-Therneau type with g(t) = 1 - KM(t-) used by
survival::cox.zph(transform = "km").

Run: python3 -I expected_km_cox.py
"""
import numpy as np
from scipy import stats

# Test data (the same vectors are written in the R tests)
CTRL_T = [2, 3, 3, 5, 6, 7, 7, 9, 11, 12, 14, 15, 18]
CTRL_E = [1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0]
DRUG_T = [3, 4, 6, 7, 8, 10, 12, 13, 16, 17, 20, 22]
DRUG_E = [1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0]
AGE = [61, 55, 70, 48, 66, 59, 73, 52, 64, 68, 57, 62, 71,
       58, 63, 49, 67, 54, 72, 60, 65, 51, 69, 56, 74]
SEX = [1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0,
       0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]
REGION = ["A", "B", "A", "A", "B", "B", "A", "B", "A", "B", "A", "B", "A",
          "B", "A", "B", "A", "B", "A", "A", "B", "A", "B", "A", "B"]

Z = stats.norm.ppf(0.975)


def km(times, events):
    t = np.asarray(times, float)
    e = np.asarray(events, int)
    ut = np.unique(t[e == 1])
    rows = []
    s = 1.0
    v = 0.0
    for u in ut:
        n = np.sum(t >= u)
        d = np.sum((t == u) & (e == 1))
        s *= 1.0 - d / n
        v += d / (n * (n - d))
        rows.append((u, n, d, s, v))
    return rows


def km_at(rows, x):
    s, v = 1.0, 0.0
    for (u, n, d, su, vu) in rows:
        if u <= x:
            s, v = su, vu
    return s, v


def ci(s, v, kind):
    se_log = np.sqrt(v)
    if s == 1.0:
        return 1.0, 1.0
    if kind == "log":
        return s * np.exp(-Z * se_log), min(1.0, s * np.exp(Z * se_log))
    if kind == "log-log":
        ll = np.log(-np.log(s))
        h = Z * se_log / abs(np.log(s))
        return np.exp(-np.exp(ll + h)), np.exp(-np.exp(ll - h))
    if kind == "plain":
        se = s * se_log
        return max(0.0, s - Z * se), min(1.0, s + Z * se)
    raise ValueError(kind)


def first_time_below(rows, which, kind, p=0.5):
    for (u, n, d, s, v) in rows:
        val = {"est": s, "lower": ci(s, v, kind)[0],
               "upper": ci(s, v, kind)[1]}[which]
        assert abs(val - p) > 1e-8, "curve equal to p: midpoint rule needed"
        if val <= p:
            return u
    return np.nan


def efron_terms(beta, X, t, e, strata):
    """Log partial likelihood, score and information (Efron ties)."""
    X = np.atleast_2d(X)
    p = X.shape[1]
    ll = 0.0
    U = np.zeros(p)
    I = np.zeros((p, p))
    eta = X @ beta
    w = np.exp(eta)
    for s in np.unique(strata):
        ins = strata == s
        for u in np.unique(t[ins & (e == 1)]):
            risk = ins & (t >= u)
            dead = ins & (t == u) & (e == 1)
            d = dead.sum()
            r0 = w[risk].sum()
            r1 = (w[risk][:, None] * X[risk]).sum(0)
            r2 = (w[risk][:, None, None] * X[risk][:, :, None] * X[risk][:, None, :]).sum(0)
            d0 = w[dead].sum()
            d1 = (w[dead][:, None] * X[dead]).sum(0)
            d2 = (w[dead][:, None, None] * X[dead][:, :, None] * X[dead][:, None, :]).sum(0)
            ll += eta[dead].sum()
            U += X[dead].sum(0)
            for k in range(d):
                f = k / d
                a0 = r0 - f * d0
                a1 = r1 - f * d1
                a2 = r2 - f * d2
                ll -= np.log(a0)
                m = a1 / a0
                U -= m
                I += a2 / a0 - np.outer(m, m)
    return ll, U, I


def cox(X, t, e, strata=None):
    X = np.atleast_2d(np.asarray(X, float))
    if X.shape[0] == 1:
        X = X.T
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    strata = np.zeros(len(t)) if strata is None else np.asarray(strata)
    beta = np.zeros(X.shape[1])
    for _ in range(100):
        ll, U, I = efron_terms(beta, X, t, e, strata)
        step = np.linalg.solve(I, U)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-13:
            break
    ll, U, I = efron_terms(beta, X, t, e, strata)
    se = np.sqrt(np.diag(np.linalg.inv(I)))
    return beta, se


def zph_test(x, t, e):
    """Score test for x * g(t), g(t) = 1 - KM(t-) centred over events."""
    x = np.asarray(x, float)
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    beta, _ = cox(x, t, e)
    rows = km(t, e)
    g_obs = np.array([1.0 - km_at(rows, ti - 1e-9)[0] for ti in t])
    g = g_obs - g_obs[e == 1].mean()
    w = np.exp(x * beta[0])
    U = 0.0
    I = np.zeros((2, 2))
    for u in np.unique(t[e == 1]):
        risk = t >= u
        dead = (t == u) & (e == 1)
        gu = g[dead][0]
        d = dead.sum()
        r0, r1, r2 = w[risk].sum(), (w * x)[risk].sum(), (w * x * x)[risk].sum()
        d0, d1, d2 = w[dead].sum(), (w * x)[dead].sum(), (w * x * x)[dead].sum()
        U += gu * x[dead].sum()
        for k in range(d):
            f = k / d
            a0, a1, a2 = r0 - f * d0, r1 - f * d1, r2 - f * d2
            m = a1 / a0
            vv = a2 / a0 - m * m
            U -= gu * m
            I += vv * np.array([[1.0, gu], [gu, gu * gu]])
    u = np.array([0.0, U])
    chisq = float(u @ np.linalg.solve(I, u))
    return beta[0], chisq, stats.chi2.sf(chisq, 1)


def fmt(x):
    return repr(float(x))


if __name__ == "__main__":
    rows_c = km(CTRL_T, CTRL_E)
    rows_d = km(DRUG_T, DRUG_E)

    print("== summarise_data: median and 95% CI (log) ==")
    for name, rows in (("Control", rows_c), ("Drug", rows_d)):
        print(name, [first_time_below(rows, w, "log") for w in ("est", "lower", "upper")])
    pooled = km(CTRL_T + DRUG_T, CTRL_E + DRUG_E)
    print("Overall", [first_time_below(pooled, w, "log") for w in ("est", "lower", "upper")])

    print("== km_table: interval 5 ==")
    for name, (tt, ee, rows) in (("Control", (CTRL_T, CTRL_E, rows_c)),
                                 ("Drug", (DRUG_T, DRUG_E, rows_d))):
        tt = np.array(tt, float)
        ee = np.array(ee, int)
        grid = np.arange(0, tt.max() + 1e-9, 5.0)
        prev = -np.inf
        for g in grid:
            s, v = km_at(rows, g)
            n_risk = int(np.sum(tt >= g))
            n_ev = int(np.sum((tt > prev) & (tt <= g) & (ee == 1)))
            n_ce = int(np.sum((tt > prev) & (tt <= g) & (ee == 0)))
            prev = g
            out = [name, g, n_risk, n_ev, n_ce, fmt(s), fmt(s * np.sqrt(v))]
            for kind in ("log", "log-log", "plain"):
                out += [fmt(c) for c in ci(s, v, kind)]
            print(out)

    t_all = np.array(CTRL_T + DRUG_T, float)
    e_all = np.array(CTRL_E + DRUG_E, int)
    arm = np.array([0] * len(CTRL_T) + [1] * len(DRUG_T), float)
    age = np.array(AGE, float)
    sex = np.array(SEX, float)
    region = np.array(REGION)

    print("== estimate_hazard_ratio (Drug vs Control) ==")
    for label, X, st in (("unadjusted, unstratified", arm[:, None], None),
                         ("unadjusted, stratified", arm[:, None], region),
                         ("adjusted, unstratified", np.column_stack([arm, age, sex]), None),
                         ("adjusted, stratified", np.column_stack([arm, age, sex]), region)):
        b, se = cox(X, t_all, e_all, st)
        hr = np.exp(b[0])
        lo, hi = np.exp(b[0] - Z * se[0]), np.exp(b[0] + Z * se[0])
        p = 2 * stats.norm.sf(abs(b[0] / se[0]))
        print(label, fmt(hr), fmt(lo), fmt(hi), fmt(p))

    print("== test_proportional_hazards (Drug vs Control, transform = km) ==")
    b, chisq, p = zph_test(arm, t_all, e_all)
    print(fmt(b), fmt(chisq), fmt(p))
