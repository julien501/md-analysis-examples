#!/usr/bin/env python3
"""Statistics for replicated MD simulations, with the replica (one run) as the unit, never the frame."""
import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 0.05
KT_KCAL = 0.001987204 * 310.0          # kT at 310 K in kcal/mol


def describe(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return {"n": len(v), "mean": float(v.mean()) if len(v) else np.nan,
            "sd": float(v.std(ddof=1)) if len(v) > 1 else np.nan}


def variant_anova(values, groups):
    """One-way ANOVA over the systems with the replicas as observations.
    icc = ICC(1), the share of the variance that comes from the system. sd_replica = pooled SD within systems."""
    df = pd.DataFrame({"v": np.asarray(values, dtype=float), "g": np.asarray(groups)}).dropna()
    g = df.groupby("g").v
    n_i, means = g.size().to_numpy(), g.mean().to_numpy()
    k, N = len(n_i), len(df)
    empty = {"icc": np.nan, "p": np.nan, "F": np.nan, "df1": k - 1, "df2": N - k, "sd_replica": np.nan,
             "n_variants": k}
    if k < 2 or N - k < 1:
        return empty
    grand = df.v.mean()
    # a constant descriptor still has rounding residuals (~1e-17), which the F-test would call highly significant
    scale = max(abs(float(grand)), float(np.abs(df.v).max()), 1.0)
    if float(df.v.max() - df.v.min()) <= 1e-12 * scale:
        return {**empty, "sd_replica": 0.0}
    ssb = float((n_i * (means - grand) ** 2).sum())
    ssw = float(((df.v - g.transform("mean")) ** 2).sum())
    msb, msw = ssb / (k - 1), ssw / (N - k)
    n0 = (N - (n_i ** 2).sum() / N) / (k - 1)          # mean group size, so unequal replica numbers are fine
    if msw > 0:
        F = msb / msw
        p = float(stats.f.sf(F, k - 1, N - k))
    else:
        F, p = (np.inf, 0.0) if msb > 0 else (np.nan, np.nan)
    denom = msb + (n0 - 1) * msw
    icc = (msb - msw) / denom if denom > 0 else np.nan
    return {"icc": float(icc), "p": p, "F": float(F), "df1": k - 1, "df2": N - k,
            "sd_replica": float(np.sqrt(msw)), "n_variants": k}


def dunnett_vs_wt(values, groups, control="WT"):
    """Each system vs the control (Dunnett's test). diff = mean difference, d = diff / pooled replica SD."""
    df = pd.DataFrame({"v": np.asarray(values, dtype=float), "g": np.asarray(groups)}).dropna()
    empty = pd.DataFrame(columns=["variant", "diff", "d", "p"])
    others = [x for x in df.g.unique() if x != control and (df.g == x).sum() >= 2]
    if control not in set(df.g) or not others:
        return empty
    ctrl = df.v[df.g == control].to_numpy()
    samples = [df.v[df.g == x].to_numpy() for x in others]
    res = stats.dunnett(*samples, control=ctrl)
    sd = variant_anova(df.v, df.g)["sd_replica"]
    return pd.DataFrame({"variant": others,
                         "diff": [s.mean() - ctrl.mean() for s in samples],
                         "d": [(s.mean() - ctrl.mean()) / sd if sd > 0 else np.nan for s in samples],
                         "p": res.pvalue})


def variant_test(values, groups, control="WT"):
    """ANOVA first. Dunnett only if the systems differ at all (p < 0.05)."""
    a = variant_anova(values, groups)
    if np.isfinite(a["p"]) and a["p"] < ALPHA:
        return a, dunnett_vs_wt(values, groups, control)
    return a, pd.DataFrame(columns=["variant", "diff", "d", "p"])


def strength(r):
    # verbal scale of Schober et al. 2018, Anesth Analg 126:1763
    if not np.isfinite(r):
        return "n/a"
    a = abs(r)
    return ("negligible" if a < 0.1 else "weak" if a < 0.4 else "moderate" if a < 0.7
            else "strong" if a < 0.9 else "very strong")


def spearman(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return {"r": np.nan, "p": np.nan, "n": int(ok.sum()), "strength": "n/a"}
    r, p = stats.spearmanr(x[ok], y[ok])
    return {"r": float(r), "p": float(p), "n": int(ok.sum()), "strength": strength(r)}


def free_energy(density):
    """-kT ln(P / P_max) in kcal/mol, 0 at the most populated value."""
    d = np.asarray(density, dtype=float)
    with np.errstate(divide="ignore"):
        return -KT_KCAL * np.log(d / d.max())


def replicas_needed(sd, delta, alpha=ALPHA, power=0.8):
    """Replicas per system to detect a difference delta (two-sided t-test, normal approximation)."""
    if not (np.isfinite(sd) and np.isfinite(delta)) or delta == 0:
        return np.nan
    z = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    return int(np.ceil(2 * z ** 2 * sd ** 2 / delta ** 2))


def bh(p):
    """Benjamini-Hochberg q values. NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(len(pv))
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q


def fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    return f"{p:.2g}" if p >= 0.001 else f"{p:.0e}"
