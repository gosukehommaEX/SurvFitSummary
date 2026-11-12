"""Independent reference values for the model fitting tests of SurvFitSummary
(tests/testthat/test-fit_survival_models.R, test-parameter_table.R,
test-covariance_table.R, test-gof_table.R, test-predict_survival.R,
test-survival_summary.R and test-spline_knot_table.R).

The computations use numpy and scipy only (no survival analysis library). The
twelve parametric distributions are written in the parameterisation of the R
package flexsurv (Jackson 2016, J Stat Softw 70(8)), with the arm effect on
the transformed location parameter; the Royston-Parmar models use the natural
cubic spline basis of Royston and Parmar (2002) on log time with the hazard,
odds and normal scales, the arm effect on gamma0, and internal knots at
equally spaced quantiles (type 7) of the log event times, the boundary knots
at their minimum and maximum. The maximum likelihood estimates are found
from several starting values (BFGS and Nelder-Mead) and refined by Newton
steps; the Hessian of minus the log-likelihood is computed by central
differences with Richardson extrapolation (initial step 1% of the
parameter, checked against half and double steps), its inverse is the
covariance matrix on the transformed scale and the delta method gives the
covariance matrix on the natural scale. Medians are found by root finding,
means and restricted means by numerical integration (means also from closed
forms where they exist), Kaplan-Meier estimates with Greenwood standard
errors.

Input: ../../extdata/fit_example.csv (written by make_fit_example.py).
Output: expected/parameters.csv, covariance.csv, gof.csv, predictions.csv,
summary.csv and knots.csv, read by the tests.

Run: python3 -I expected_fits.py   (about two minutes; prints the check of
the numerical Hessian)
"""
import numpy as np
from scipy import special, optimize, integrate, stats

Z_975 = stats.norm.ppf(0.975)

# ---------------------------------------------------------------------------
# Distributions: names of the parameters, which are log-transformed, which
# one is the location parameter (where the arm effect enters, on the
# transformed scale), and log density / log survival on the natural scale.
# ---------------------------------------------------------------------------


def _gengamma_logf_logs(t, mu, sigma, q):
    lt = np.log(t)
    w = (lt - mu) / sigma
    q = np.broadcast_to(q, np.shape(w)).astype(float)
    small = np.abs(q) < 1e-8
    out_f = np.empty_like(w)
    out_s = np.empty_like(w)
    if np.any(small):
        ws = w[small]
        sg = np.broadcast_to(sigma, np.shape(w))[small]
        out_f[small] = stats.norm.logpdf(ws) - np.log(sg) - lt[small]
        out_s[small] = stats.norm.logsf(ws)
    big = ~small
    if np.any(big):
        qb = q[big]
        wb = w[big]
        sg = np.broadcast_to(sigma, np.shape(w))[big]
        g = qb ** -2
        qw = qb * wb
        out_f[big] = (-np.log(sg * t[big]) + np.log(np.abs(qb)) + g * np.log(g)
                      + g * (qw - np.exp(qw)) - special.gammaln(g))
        u = g * np.exp(qw)
        s = np.where(qb > 0, special.gammaincc(g, u), special.gammainc(g, u))
        out_s[big] = np.log(s)
    return out_f, out_s


def _genf_orig_logf_logs(t, mu, sigma, s1, s2):
    lt = np.log(t)
    w = (lt - mu) / sigma
    logf = (s1 * (np.log(s1) - np.log(s2)) + s1 * w - np.log(sigma) - lt
            - (s1 + s2) * np.log1p(s1 * np.exp(w) / s2) - special.betaln(s1, s2))
    x = s2 / (s2 + s1 * np.exp(w))
    logs = np.log(special.betainc(s2, s1, x))
    return logf, logs


def _genf_logf_logs(t, mu, sigma, q, p):
    tmp = q ** 2 + 2 * p
    delta = np.sqrt(tmp)
    s1 = 2 / (tmp + q * delta)
    s2 = 2 / (tmp - q * delta)
    return _genf_orig_logf_logs(t, mu, sigma / delta, s1, s2)


def dist_exp(t, rate):
    return np.log(rate) - rate * t, -rate * t


def dist_weibull(t, shape, scale):
    z = (t / scale) ** shape
    logf = np.log(shape) - np.log(scale) + (shape - 1) * (np.log(t) - np.log(scale)) - z
    return logf, -z


def dist_weibullph(t, shape, scale):
    z = scale * t ** shape
    logf = np.log(shape) + np.log(scale) + (shape - 1) * np.log(t) - z
    return logf, -z


def dist_lnorm(t, meanlog, sdlog):
    w = (np.log(t) - meanlog) / sdlog
    return stats.norm.logpdf(w) - np.log(sdlog) - np.log(t), stats.norm.logsf(w)


def dist_llogis(t, shape, scale):
    u = (t / scale) ** shape
    logf = (np.log(shape) - np.log(scale) + (shape - 1) * (np.log(t) - np.log(scale))
            - 2 * np.log1p(u))
    return logf, -np.log1p(u)


def dist_gompertz(t, shape, rate):
    shape = np.broadcast_to(shape, np.shape(t)).astype(float)
    at = shape * t
    cum = np.where(np.abs(shape) < 1e-12, rate * t,
                   rate * np.expm1(at) / np.where(np.abs(shape) < 1e-12, 1.0, shape))
    logh = np.log(rate) + at
    return logh - cum, -cum


def dist_gamma(t, shape, rate):
    logf = shape * np.log(rate) + (shape - 1) * np.log(t) - rate * t - special.gammaln(shape)
    return logf, np.log(special.gammaincc(shape, rate * t))


def dist_gengamma(t, mu, sigma, q):
    return _gengamma_logf_logs(np.asarray(t, float), mu, sigma, q)


def dist_gengamma_orig(t, shape, scale, k):
    u = (t / scale) ** shape
    logf = (np.log(shape) + (shape * k - 1) * np.log(t) - shape * k * np.log(scale)
            - special.gammaln(k) - u)
    return logf, np.log(special.gammaincc(k, u))


def dist_genf(t, mu, sigma, q, p):
    return _genf_logf_logs(np.asarray(t, float), mu, sigma, q, p)


def dist_genf_orig(t, mu, sigma, s1, s2):
    return _genf_orig_logf_logs(np.asarray(t, float), mu, sigma, s1, s2)


def dist_rayleigh(t, scale):
    return dist_weibull(t, 2.0, scale)


# name: (function, parameter names, log-transformed flags, location index)
DISTS = {
    "exp": (dist_exp, ["rate"], [True], 0),
    "weibull": (dist_weibull, ["shape", "scale"], [True, True], 1),
    "weibullPH": (dist_weibullph, ["shape", "scale"], [True, True], 1),
    "rayleigh": (dist_rayleigh, ["scale"], [True], 0),
    "lnorm": (dist_lnorm, ["meanlog", "sdlog"], [False, True], 0),
    "llogis": (dist_llogis, ["shape", "scale"], [True, True], 1),
    "gompertz": (dist_gompertz, ["shape", "rate"], [False, True], 1),
    "gamma": (dist_gamma, ["shape", "rate"], [True, True], 1),
    "gengamma": (dist_gengamma, ["mu", "sigma", "Q"], [False, True, False], 0),
    "gengamma.orig": (dist_gengamma_orig, ["shape", "scale", "k"], [True, True, True], 1),
    "genf": (dist_genf, ["mu", "sigma", "Q", "P"], [False, True, False, True], 0),
    "genf.orig": (dist_genf_orig, ["mu", "sigma", "s1", "s2"], [False, True, True, True], 0),
}


# ---------------------------------------------------------------------------
# Royston-Parmar natural cubic spline basis (flexsurv spline = "rp")
# ---------------------------------------------------------------------------


def rp_basis(knots, x):
    """Basis [1, x, v_1(x), ..., v_m(x)] and its derivative with respect to x."""
    x = np.asarray(x, float)
    kmin, kmax = knots[0], knots[-1]
    internal = knots[1:-1]
    b = [np.ones_like(x), x.copy()]
    db = [np.zeros_like(x), np.ones_like(x)]
    for kj in internal:
        lam = (kmax - kj) / (kmax - kmin)
        b.append(np.maximum(x - kj, 0) ** 3 - lam * np.maximum(x - kmin, 0) ** 3
                 - (1 - lam) * np.maximum(x - kmax, 0) ** 3)
        db.append(3 * np.maximum(x - kj, 0) ** 2 - 3 * lam * np.maximum(x - kmin, 0) ** 2
                  - 3 * (1 - lam) * np.maximum(x - kmax, 0) ** 2)
    return np.column_stack(b), np.column_stack(db)


def spline_logf_logs(t, gamma, knots, scale, shift):
    t = np.asarray(t, float)
    x = np.log(t)
    b, db = rp_basis(knots, x)
    s = b @ gamma + shift
    dsdx = db @ gamma
    with np.errstate(invalid="ignore", divide="ignore"):
        logd = np.log(dsdx) - x
    if scale == "hazard":
        logs = -np.exp(s)
        logf = s - np.exp(s) + logd
    elif scale == "odds":
        logs = -np.logaddexp(0, s)
        logf = s - 2 * np.logaddexp(0, s) + logd
    elif scale == "normal":
        logs = stats.norm.logcdf(-s)
        logf = stats.norm.logpdf(s) + logd
    else:
        raise ValueError(scale)
    return logf, logs


def default_knots(time, event, k):
    lt = np.log(np.asarray(time, float)[np.asarray(event) == 1])
    probs = np.linspace(0, 1, k + 2)[1:-1]
    internal = np.quantile(lt, probs)  # type 7, the default of R quantile()
    return np.concatenate([[lt.min()], internal, [lt.max()]]), probs


# ---------------------------------------------------------------------------
# Model object: theta is on the transformed (real-line) scale, ordered as in
# flexsurv: baseline parameters, then the arm effects.
# ---------------------------------------------------------------------------


class Model:
    def __init__(self, family, time, event, arm_index=None, n_arms=1,
                 dist=None, scale=None, knots=None):
        self.family = family
        self.time = np.asarray(time, float)
        self.event = np.asarray(event, int)
        self.n_cov = n_arms - 1
        if arm_index is None:
            arm_index = np.zeros(len(self.time), int)
        self.z = np.zeros((len(self.time), self.n_cov))
        for j in range(self.n_cov):
            self.z[:, j] = (np.asarray(arm_index) == j + 1).astype(float)
        if family == "parametric":
            self.dist = dist
            f, names, logflag, loc = DISTS[dist]
            self.fun, self.base_names, self.logflag, self.loc = f, names, logflag, loc
        else:
            self.dist = "spline"
            self.scale = scale
            self.knots = np.asarray(knots, float)
            m = len(self.knots) - 2
            self.base_names = [f"gamma{i}" for i in range(m + 2)]
            self.logflag = [False] * (m + 2)
            self.loc = 0
        self.npars = len(self.base_names) + self.n_cov

    def natural_base(self, theta):
        nb = len(self.base_names)
        return np.array([np.exp(v) if lg else v
                         for v, lg in zip(theta[:nb], self.logflag)])

    def _pars_per_obs(self, theta, z):
        nb = len(self.base_names)
        base_t = np.asarray(theta[:nb], float)
        beta = np.asarray(theta[nb:], float)
        shift = z @ beta if self.n_cov > 0 else np.zeros(z.shape[0])
        pars = []
        for i in range(nb):
            v = base_t[i] + (shift if i == self.loc else 0.0)
            pars.append(np.exp(v) if self.logflag[i] else v)
        return pars, shift

    def logf_logs(self, theta, t, z):
        if self.family == "parametric":
            pars, _ = self._pars_per_obs(theta, z)
            return self.fun(t, *pars)
        nb = len(self.base_names)
        gamma = np.asarray(theta[:nb], float)
        shift = z @ np.asarray(theta[nb:], float) if self.n_cov > 0 else 0.0
        return spline_logf_logs(t, gamma, self.knots, self.scale, shift)

    def loglik(self, theta):
        with np.errstate(all="ignore"):
            logf, logs = self.logf_logs(theta, self.time, self.z)
            ll = np.sum(np.where(self.event == 1, logf, logs))
        return ll if np.isfinite(ll) else -np.inf

    def nll(self, theta):
        ll = self.loglik(theta)
        return -ll if np.isfinite(ll) else 1e300

    # predictions for one arm (row of z)
    def surv(self, theta, t, z_row):
        t = np.atleast_1d(np.asarray(t, float))
        zz = np.tile(np.asarray(z_row, float), (len(t), 1))
        with np.errstate(all="ignore"):
            _, logs = self.logf_logs(theta, t, zz)
        return np.exp(logs)

    def hazard(self, theta, t, z_row):
        t = np.atleast_1d(np.asarray(t, float))
        zz = np.tile(np.asarray(z_row, float), (len(t), 1))
        logf, logs = self.logf_logs(theta, t, zz)
        return np.exp(logf - logs)


# ---------------------------------------------------------------------------
# Numerical derivatives with Richardson extrapolation (as in numDeriv). The
# initial step is 1% of the parameter (at least 1e-4) and is halved 3 times.
# Much smaller steps amplify the rounding noise of the log-likelihoods with
# incomplete beta functions (generalized F: up to 6% error in the covariance
# matrix with an initial step of 0.1%), and much larger steps leave the region
# where the derivative of a spline is positive. main() checks every covariance
# matrix against those with half and double the initial step.
# ---------------------------------------------------------------------------


def gradient(f, x, d=1e-2, r=4, v=2):
    x = np.asarray(x, float)
    g = np.zeros_like(x)
    for i in range(len(x)):
        h0 = d * max(abs(x[i]), 1e-2)
        a = []
        for k in range(r):
            h = h0 / v ** k
            e = np.zeros_like(x)
            e[i] = h
            a.append((f(x + e) - f(x - e)) / (2 * h))
        a = np.array(a)
        for m in range(1, r):
            a = (a[1:] * 4 ** m - a[:-1]) / (4 ** m - 1)
        g[i] = a[0]
    return g


def hessian(f, x, d=1e-2, r=4, v=2):
    x = np.asarray(x, float)
    n = len(x)
    hmat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            hi0 = d * max(abs(x[i]), 1e-2)
            hj0 = d * max(abs(x[j]), 1e-2)
            a = []
            for k in range(r):
                hi = hi0 / v ** k
                hj = hj0 / v ** k
                ei = np.zeros(n)
                ej = np.zeros(n)
                ei[i] = hi
                ej[j] = hj
                val = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * hi * hj)
                a.append(val)
            a = np.array(a)
            for m in range(1, r):
                a = (a[1:] * 4 ** m - a[:-1]) / (4 ** m - 1)
            hmat[i, j] = hmat[j, i] = a[0]
    return hmat


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------


def fit(model, starts, polish=True):
    """Maximum likelihood: BFGS, Nelder-Mead and BFGS again from each start,
    the best result refined by Newton steps with a numerical Hessian.
    Returns the estimate, the covariance matrix (inverse Hessian of minus the
    log-likelihood), the log-likelihood, the gradient and the Hessian."""
    best = None
    for s in starts:
        s = np.asarray(s, float)
        if not np.isfinite(model.loglik(s)):
            continue
        res = optimize.minimize(model.nll, s, method="BFGS",
                                options={"gtol": 1e-8, "maxiter": 5000})
        res = optimize.minimize(model.nll, res.x, method="Nelder-Mead",
                                options={"xatol": 1e-9, "fatol": 1e-11,
                                         "maxiter": 40000, "maxfev": 80000})
        res = optimize.minimize(model.nll, res.x, method="BFGS",
                                options={"gtol": 1e-9, "maxiter": 5000})
        if best is None or res.fun < best.fun - 1e-9:
            best = res
    theta = best.x.copy()
    if polish:
        for _ in range(30):
            g = gradient(model.nll, theta)
            hm = hessian(model.nll, theta)
            if not np.all(np.linalg.eigvalsh(hm) > 0):
                break
            step = np.linalg.solve(hm, g)
            # A Newton step is accepted up to the rounding noise of the
            # log-likelihood; otherwise it is halved.
            lam = 1.0
            while lam > 1e-4 and model.nll(theta - lam * step) > model.nll(theta) + 1e-9:
                lam /= 2
            if lam <= 1e-4:
                break
            theta = theta - lam * step
            if np.max(np.abs(lam * step)) < 1e-10:
                break
    hm = hessian(model.nll, theta)
    cov = np.linalg.inv(hm)
    return theta, cov, model.loglik(theta), gradient(model.nll, theta), hm


def default_starts(model, rng=None, n_random=6):
    """Starting values on the transformed scale (several, to find the global
    maximum)."""
    t = model.time
    e = model.event
    lt = np.log(t)
    m, s = np.mean(lt), np.std(lt)
    rate = e.sum() / t.sum()
    base = {
        "exp": [np.log(rate)],
        "weibull": [0.0, m],
        "weibullPH": [0.0, np.log(rate)],
        "rayleigh": [np.log(np.sqrt(np.sum(t ** 2) / e.sum()))],
        "lnorm": [m, np.log(s)],
        "llogis": [0.0, m],
        "gompertz": [0.0, np.log(rate)],
        "gamma": [0.0, np.log(rate)],
        "gengamma": [m, np.log(s), 0.5],
        "gengamma.orig": [0.0, m, 0.0],
        "genf": [m, np.log(s), 0.5, 0.0],
        "genf.orig": [m, np.log(s), 0.0, 0.0],
    }
    if model.family == "parametric":
        b0 = np.array(base[model.dist], float)
    else:
        b0 = np.zeros(len(model.base_names))
        b0[0] = np.log(-np.log(0.5)) - m if model.scale == "hazard" else -m
        b0[1] = 1.0
    starts = [np.concatenate([b0, np.zeros(model.n_cov)])]
    rng = rng or np.random.default_rng(2026)
    for _ in range(n_random):
        starts.append(starts[0] + rng.normal(0, 0.3, len(starts[0])))
    return starts


def natural_summary(model, theta, cov):
    """Estimates, confidence limits and standard errors on both scales, and
    the covariance matrix on the natural scale (delta method)."""
    nb = len(model.base_names)
    logflag = np.array(list(model.logflag) + [False] * model.n_cov)
    se_t = np.sqrt(np.diag(cov))
    lcl_t = theta - Z_975 * se_t
    ucl_t = theta + Z_975 * se_t
    est_n = np.where(logflag, np.exp(theta), theta)
    lcl_n = np.where(logflag, np.exp(lcl_t), lcl_t)
    ucl_n = np.where(logflag, np.exp(ucl_t), ucl_t)
    jac = np.where(logflag, est_n, 1.0)
    cov_n = cov * np.outer(jac, jac)
    se_n = np.sqrt(np.diag(cov_n))
    return {"est_t": theta, "lcl_t": lcl_t, "ucl_t": ucl_t, "se_t": se_t,
            "est_n": est_n, "lcl_n": lcl_n, "ucl_n": ucl_n, "se_n": se_n,
            "cov_t": cov, "cov_n": cov_n,
            "chol_t": np.linalg.cholesky(cov), "chol_n": np.linalg.cholesky(cov_n)}


def median_time(model, theta, z_row):
    f = lambda lt: model.surv(theta, np.exp(lt), z_row)[0] - 0.5
    return float(np.exp(optimize.brentq(f, -20, 20, xtol=1e-14, rtol=1e-14)))


def mean_time(model, theta, z_row):
    f = lambda t: model.surv(theta, t, z_row)[0]
    med = median_time(model, theta, z_row)
    pieces = [0, med, 10 * med, 100 * med]
    total = 0.0
    for a, b in zip(pieces[:-1], pieces[1:]):
        total += integrate.quad(f, a, b, epsabs=1e-13, epsrel=1e-12, limit=500)[0]
    tail = integrate.quad(f, pieces[-1], np.inf, epsabs=1e-13, epsrel=1e-12, limit=500)[0]
    return total + tail, f(100 * med)


def rmst(model, theta, z_row, tau):
    f = lambda t: model.surv(theta, t, z_row)[0]
    return integrate.quad(f, 0, tau, epsabs=1e-13, epsrel=1e-12, limit=500)[0]


def km(time, event):
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    rows = []
    s, v = 1.0, 0.0
    for u in np.unique(t[e == 1]):
        n = np.sum(t >= u)
        d = np.sum((t == u) & (e == 1))
        s *= 1.0 - d / n
        v += d / (n * (n - d)) if n > d else np.inf
        rows.append((u, s, v))
    return rows


def km_at(time, event, x):
    s, v = 1.0, 0.0
    for (u, su, vu) in km(time, event):
        if u <= x:
            s, v = su, vu
    return s, s * np.sqrt(v)


# ---------------------------------------------------------------------------
# Fitting all the distributions of one group with related starting values
# ---------------------------------------------------------------------------

ORDER = ["exp", "weibull", "weibullPH", "rayleigh", "lnorm", "llogis", "gompertz",
         "gamma", "gengamma", "gengamma.orig", "genf", "genf.orig"]


def related_starts(dist, fits, n_cov):
    """Starting values built from the fits of simpler distributions."""
    out = []
    zeros = np.zeros(n_cov)
    gg = fits.get("gengamma")
    gf = fits.get("genf")
    if dist == "gengamma.orig" and gg is not None:
        mu, ls, q = gg[0][:3]
        if q > 0:
            sigma = np.exp(ls)
            k = q ** -2
            b = q / sigma
            out.append(np.r_[np.log(b), mu - np.log(k) / b, np.log(k), gg[0][3:]])
    if dist == "genf" and gg is not None:
        for lp in [-3.0, -1.0, 0.0, 1.0]:
            out.append(np.r_[gg[0][:3], lp, gg[0][3:]])
    if dist == "genf.orig" and gf is not None:
        mu, ls, q, lp = gf[0][:4]
        p = np.exp(lp)
        tmp = q ** 2 + 2 * p
        delta = np.sqrt(tmp)
        s1 = 2 / (tmp + q * delta)
        s2 = 2 / (tmp - q * delta)
        out.append(np.r_[mu, ls - np.log(delta), np.log(s1), np.log(s2), gf[0][4:]])
    return out


def fit_group(time, event, arm_index=None, n_arms=1, dists=ORDER, pooled_fits=None):
    fits = {}
    for d in dists:
        mod = Model("parametric", time, event, arm_index, n_arms, dist=d)
        starts = default_starts(mod) + related_starts(d, fits, mod.n_cov)
        if pooled_fits is not None and d in pooled_fits:
            starts.append(np.r_[pooled_fits[d][0], np.zeros(mod.n_cov)])
        theta, cov, ll, g, hm = fit(mod, starts)
        fits[d] = (theta, cov, ll, g, hm, mod)
    return fits


SPLINE_SCALES = ["hazard", "odds", "normal"]


def spline_starts(scale, n_internal, fits, n_cov):
    """Starting values from the k = 0 special cases: Weibull (hazard),
    log-logistic (odds) and log-normal (normal)."""
    if scale == "hazard":
        a, m = np.exp(fits["weibullPH"][0][:2])
        g0, g1 = np.log(m), a
    elif scale == "odds":
        a, b = np.exp(fits["llogis"][0][:2])
        g0, g1 = -a * np.log(b), a
    else:
        mu, s = fits["lnorm"][0][0], np.exp(fits["lnorm"][0][1])
        g0, g1 = -mu / s, 1 / s
    base = np.r_[g0, g1, np.zeros(n_internal)]
    starts = [np.r_[base, np.zeros(n_cov)]]
    rng = np.random.default_rng(7)
    for _ in range(4):
        starts.append(np.r_[base + rng.normal(0, [0.1, 0.05] + [0.005] * n_internal),
                            np.zeros(n_cov)])
    return starts


def fit_splines(time, event, arm_index=None, n_arms=1, n_knots=(1, 2, 3), fits=None):
    out = {}
    for scale in SPLINE_SCALES:
        for k in n_knots:
            knots, probs = default_knots(time, event, k)
            mod = Model("spline", time, event, arm_index, n_arms, scale=scale, knots=knots)
            starts = spline_starts(scale, k, fits, mod.n_cov)
            theta, cov, ll, g, hm = fit(mod, starts)
            out[(scale, k)] = (theta, cov, ll, g, hm, mod, probs)
    return out


# ---------------------------------------------------------------------------
# Expected values for the tests
# ---------------------------------------------------------------------------

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "extdata", "fit_example.csv")
OUT = os.path.join(HERE, "expected")

LABELS = {"exp": "Exponential", "weibull": "Weibull", "weibullPH": "Weibull (PH)",
          "rayleigh": "Rayleigh", "lnorm": "Lognormal", "llogis": "Log-logistic",
          "gompertz": "Gompertz", "gamma": "Gamma", "gengamma": "Generalized Gamma",
          "gengamma.orig": "Generalized Gamma (original)", "genf": "Generalized F",
          "genf.orig": "Generalized F (original)"}
SCALE_LABELS = {"hazard": "Hazards", "odds": "Odds", "normal": "Normal"}
TIMES = [0.0, 6.0, 12.0, 24.0, 36.0, 60.0, 120.0]
HORIZON = 60.0
ARM_LEVELS = ["Control", "Treatment"]


def spline_label(scale, k):
    return f"{SCALE_LABELS[scale]}, {k} knot" + ("s" if k > 1 else "")


def fmt(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return "%.12g" % float(x)


def read_data(path=DATA):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    arm = np.array([r["TRT01P"] for r in rows])
    time = np.array([float(r["AVAL"]) for r in rows])
    event = 1 - np.array([int(r["CNSR"]) for r in rows])
    return arm, time, event


def closed_form_mean(dist, nat):
    """Means with a closed form, used to check the numerical integral."""
    p = nat
    if dist == "exp":
        return 1 / p[0]
    if dist == "weibull":
        return p[1] * special.gamma(1 + 1 / p[0])
    if dist == "weibullPH":
        return p[1] ** (-1 / p[0]) * special.gamma(1 + 1 / p[0])
    if dist == "rayleigh":
        return p[0] * special.gamma(1.5)
    if dist == "lnorm":
        return np.exp(p[0] + p[1] ** 2 / 2)
    if dist == "llogis":
        return p[1] * (np.pi / p[0]) / np.sin(np.pi / p[0]) if p[0] > 1 else np.inf
    if dist == "gamma":
        return p[0] / p[1]
    if dist == "gengamma.orig":
        b, a, k = p
        return a * np.exp(special.gammaln(k + 1 / b) - special.gammaln(k))
    if dist == "gengamma":
        mu, sigma, q = p
        k = q ** -2
        b = q / sigma
        a = np.exp(mu - np.log(k) / b)
        return a * np.exp(special.gammaln(k + 1 / b) - special.gammaln(k))
    if dist in ("genf", "genf.orig"):
        if dist == "genf":
            mu, sigma, q, pp = p
            tmp = q ** 2 + 2 * pp
            delta = np.sqrt(tmp)
            s1, s2 = 2 / (tmp + q * delta), 2 / (tmp - q * delta)
            sigma = sigma / delta
        else:
            mu, sigma, s1, s2 = p
        if s2 <= sigma:
            return np.inf
        return (np.exp(mu) * (s2 / s1) ** sigma
                * np.exp(special.betaln(s1 + sigma, s2 - sigma) - special.betaln(s1, s2)))
    return None


def main():
    arm, time, event = read_data()
    os.makedirs(OUT, exist_ok=True)
    groups = [("independent", "Control", arm == "Control"),
              ("independent", "Treatment", arm == "Treatment"),
              ("dependent", "Both arms", np.ones(len(time), bool))]
    par_rows, cov_rows, gof_rows, pred_rows, sum_rows, knot_rows = [], [], [], [], [], []
    hessian_check = []
    for fit_type, arm_label, m in groups:
        t, e = time[m], event[m]
        if fit_type == "dependent":
            ai = (arm[m] == "Treatment").astype(int)
            n_arms = 2
            pooled = fit_group(t, e)
            fits = fit_group(t, e, ai, n_arms, pooled_fits=pooled)
            pred_arms = [(a, [float(j)]) for j, a in enumerate(ARM_LEVELS)]
        else:
            ai, n_arms = None, 1
            fits = fit_group(t, e)
            pred_arms = [(arm_label, [])]
        splines = fit_splines(t, e, ai, n_arms, fits=fits)
        models = [(LABELS[d], fits[d][:6], d) for d in ORDER]
        models += [(spline_label(sc, k), v[:6], "spline") for (sc, k), v in splines.items()]
        for label, (theta, cov, ll, g, hm, mod), dist in models:
            names = list(mod.base_names) + [f"ARM{a}" for a in ARM_LEVELS[1:n_arms]]
            ns = natural_summary(mod, theta, cov)
            # Stability of the covariance matrix: agreement with half or double
            # the initial step of the numerical Hessian
            scale_ij = np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
            diffs = [np.max(np.abs(np.linalg.inv(hessian(mod.nll, theta, d=d)) - cov) / scale_ij)
                     for d in (5e-3, 2e-2)]
            hessian_check.append((min(diffs), fit_type, arm_label, label))
            for scale, suf in [("natural", "n"), ("transformed", "t")]:
                for i, nm in enumerate(names):
                    par_rows.append([fit_type, arm_label, label, nm, scale,
                                     ns["est_" + suf][i], ns["lcl_" + suf][i],
                                     ns["ucl_" + suf][i], ns["se_" + suf][i]])
                # The Cholesky factor is determined by the covariance matrix, so
                # the tests check L L' = covariance instead of storing L.
                mat = ns["cov_" + suf]
                for i, ri in enumerate(names):
                    for j, cj in enumerate(names):
                        cov_rows.append([fit_type, arm_label, label, scale, ri, cj, mat[i, j]])
            npars = len(theta)
            n = len(t)
            gof_rows.append([fit_type, arm_label, label, n, int(e.sum()), npars, ll,
                             -2 * ll + 2 * npars, -2 * ll + np.log(n) * npars])
            for pa, z_row in pred_arms:
                am = arm == pa
                tmax = time[am].max()
                for tt in TIMES:
                    s = mod.surv(theta, max(tt, 0.0), z_row)[0] if tt > 0 else 1.0
                    h = mod.hazard(theta, tt, z_row)[0] if tt > 0 else None
                    k_est = km_at(time[am], event[am], tt)[0] if tt <= tmax else None
                    pred_rows.append([fit_type, pa, label, tt, s, h, k_est])
                med = median_time(mod, theta, z_row)
                mean_num, tail = mean_time(mod, theta, z_row)
                if dist != "spline":
                    nat = mod.natural_base(theta.copy())
                    if fit_type == "dependent" and z_row[0] == 1.0:
                        th2 = theta.copy()
                        th2[mod.loc] += theta[len(mod.base_names)]
                        nat = mod.natural_base(th2)
                    cf = closed_form_mean(dist, nat)
                    check = cf is not None and np.isfinite(cf) and abs(mean_num / cf - 1) < 1e-8
                    if cf is None:
                        check = tail < 1e-10
                else:
                    check = tail < 1e-10
                sum_rows.append([fit_type, pa, label, med, mean_num, rmst(mod, theta, z_row, HORIZON),
                                 int(check)])
        for k in (1, 2, 3):
            knots, probs = default_knots(t, e, k)
            for j in range(k):
                # Comparisons on the log scale: a knot at an observed log time
                # (a quantile without interpolation) is that log time exactly.
                lt = knots[j + 1]
                s, se = km_at(np.log(t), e, lt)
                knot_rows.append([fit_type, arm_label, k, j + 1, probs[j], lt, np.exp(lt),
                                  int(np.sum(np.log(t) >= lt)), s, se])
    write_csv("parameters.csv", ["fit_type", "arm", "model", "parameter", "scale", "est", "lcl",
                                 "ucl", "se"], par_rows)
    write_csv("covariance.csv", ["fit_type", "arm", "model", "scale", "row", "col", "value"],
              cov_rows)
    write_csv("gof.csv", ["fit_type", "arm", "model", "n", "events", "npars", "loglik", "AIC",
                          "BIC"], gof_rows)
    write_csv("predictions.csv", ["fit_type", "arm", "model", "time", "surv", "hazard", "km"],
              pred_rows)
    write_csv("summary.csv", ["fit_type", "arm", "model", "median", "mean", "rmst",
                              "mean_check"], sum_rows)
    write_csv("knots.csv", ["fit_type", "arm", "n_knots", "knot", "quantile", "log_time", "time",
                            "n_risk", "surv_est", "surv_se"], knot_rows)
    worst = max(hessian_check)
    print(f"Covariance matrices of {len(hessian_check)} models: largest difference from the "
          f"half or double step {worst[0]:.1e} of sqrt(V_ii V_jj) ({', '.join(worst[1:])})")
    if worst[0] > 1e-3:
        raise SystemExit("The numerical Hessian is not stable enough.")


def write_csv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\r\n")
        w.writerow(header)
        for r in rows:
            w.writerow([v if isinstance(v, str) else fmt(v) for v in r])


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()
