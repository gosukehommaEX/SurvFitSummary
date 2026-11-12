"""Simulates inst/extdata/fit_example.csv, the dataset used by the tests of
the model fitting functions of SurvFitSummary and by their examples.

Event times follow the generalized F distribution of Prentice (1975) in the
parameterisation of flexsurv (mu, sigma, Q, P), with a different location mu
in each arm: log T = mu + (sigma / delta) * log X, X ~ F(2 s1, 2 s2), where
delta = (Q^2 + 2 P)^(1/2), s1 = 2 / (Q^2 + 2 P + Q delta) and
s2 = 2 / (Q^2 + 2 P - Q delta). With Q > 0 and P > 0 the maximum likelihood
estimates of all the distributions of the package, including the generalized
gamma and generalized F in their original parameterisations, lie inside the
parameter space. Patients enter uniformly over 24 months, the data cut-off is
at 72 months and dropout is exponential with rate 0.01 per month. Times are
in months, rounded to 4 decimals. The seed was chosen among the first four so
that the generalized F has an interior maximum in both arms.

Run: python3 -I make_fit_example.py [output.csv]
"""
import sys
import numpy as np

SEED = 4
N = 300
ARMS = [("Control", 2.4), ("Treatment", 2.8)]
SIGMA, Q, P = 0.6, 0.6, 1.5
ACCRUAL, CUTOFF, DROPOUT = 24.0, 72.0, 0.01


def genf_to_orig(sigma, q, p):
    tmp = q ** 2 + 2 * p
    delta = np.sqrt(tmp)
    return sigma / delta, 2 / (tmp + q * delta), 2 / (tmp - q * delta)


def simulate(seed=SEED):
    sigma_orig, s1, s2 = genf_to_orig(SIGMA, Q, P)
    rng = np.random.default_rng(seed)
    rows = []
    k = 0
    for arm, mu in ARMS:
        x = (rng.gamma(s1, 1.0, N) / s1) / (rng.gamma(s2, 1.0, N) / s2)
        event_time = np.exp(mu + sigma_orig * np.log(x))
        entry = rng.uniform(0, ACCRUAL, N)
        dropout = rng.exponential(1 / DROPOUT, N)
        censor_time = np.minimum(CUTOFF - entry, dropout)
        time = np.round(np.minimum(event_time, censor_time), 4)
        event = (event_time <= censor_time).astype(int)
        for i in range(N):
            k += 1
            rows.append((f"P{k:03d}", arm, time[i], 1 - event[i]))
    return rows


def write(rows, path):
    with open(path, "w", newline="") as f:
        f.write("USUBJID,TRT01P,AVAL,CNSR\r\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]}\r\n")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "fit_example.csv"
    write(simulate(), out)
