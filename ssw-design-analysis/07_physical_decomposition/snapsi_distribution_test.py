#!/usr/bin/env python3
"""
snapsi_distribution_test.py
===========================
Does a sudden stratospheric warming SHIFT the surface distribution, or SPLIT it?

WHY THIS IS THE TEST THAT MATTERS
  The project's central physical claim -- "downward propagation" is a threshold
  on ONE shifted population, not a class of event -- currently rests on CMIP6
  (n=1888) under an assumed structure, and is only "consistent with" 43 observed
  events that have essentially no power. That is the weakest joint in the
  argument and the obvious line of attack.

  SNAPSI closes it EXPERIMENTALLY:
      control  no SSW forcing at all          -> the unforced distribution
      nudged   stratosphere nudged to the SAME observed SSW for every member
                                              -> the forced distribution
  Members differ only in tropospheric initial condition and internal noise, so
  the forcing is identical by construction. Comparing the two distributions
  measures what an SSW does to the SHAPE of the surface response, causally.

  one shifted population  ->  mean moves, variance and shape do not
  two classes             ->  variance inflates, or a second mode appears

WHAT IS MEASURED, per (centre, initialisation), then pooled
  shift            mean(nudged) - mean(control), in control sigma
  variance ratio   var(nudged) / var(control)   -- the causal forced variance
  pure translation KS between the two after removing each mean
  mixture          Gaussian-mixture BIC k=2 vs k=1, CALIBRATED against the
                   control (one population by construction), never against a
                   theoretical null
  power            what component separation this design could have detected

GUARDS DECIDED BEFORE RUNNING
  - standardise per (centre, init) by the CONTROL spread; the raw spread varies
    3-4x across cases and pooling raw values would be meaningless
  - the s20190108 initialisation is reported separately, at BOTH centres.
    It initialises 2019-01-08 with onset 2019-01-02, so onset PRECEDES the
    initialisation by 6 days and the +8..+25 post-onset window falls at lead
    2-19 days, before the ensemble has diverged. Its control spread is
    137-138 Pa against 205-466 Pa elsewhere -- that is forecast spread growth,
    not climatological spread, and standardising by it inflates both the shift
    and the variance ratio. Flagged in advance from the sd, then confirmed by
    the lead arithmetic; it is a property of the INITIALISATION, not of a
    centre, which is why both centres' s20190108 cases are separated.
  - NRL is excluded: nudged-minus-control ensemble effect is identically zero at
    all four initialisations and one submission duplicates control exactly.

Output: results/current/8_experiment/snapsi_distribution_test.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "8_experiment"
RESULTS.mkdir(parents=True, exist_ok=True)
RED = HERE.parents[0] / "03_data_ingestion" / "_snapsi_reduced"

ONSET = {"s20180125": "2018-02-12", "s20180208": "2018-02-12",
         "s20181213": "2019-01-02", "s20190108": "2019-01-02"}
INIT_DATE = {"s20180125": "2018-01-25", "s20180208": "2018-02-08",
             "s20181213": "2018-12-13", "s20190108": "2019-01-08"}
WINDOW = (8, 25)
SHORT_LEAD_INIT = "s20190108"   # onset PRECEDES init; see the guard note
SEED = 20260804
N_BOOT = 4000
N_POWER = 200


def load(centre, exp, init):
    fs = sorted(RED.glob(f"{centre}_{exp}_{init}_*.parquet"))
    return pd.concat([pd.read_parquet(f) for f in fs]) if fs else None


def window_means(df, init):
    off = (pd.Timestamp(ONSET[init]) - pd.Timestamp(INIT_DATE[init])).days
    w = df[(df["lead_days"] >= WINDOW[0] + off)
           & (df["lead_days"] <= WINDOW[1] + off)]
    return w.groupby("member")["psl_cap"].mean().values


def bic_k2_minus_k1(x, rng):
    """BIC(k=1) - BIC(k=2). Positive favours two components."""
    from sklearn.mixture import GaussianMixture
    x = np.asarray(x, float).reshape(-1, 1)
    if len(x) < 20:
        return np.nan
    b = []
    for k in (1, 2):
        g = GaussianMixture(k, n_init=5, random_state=int(rng.integers(1 << 31)))
        g.fit(x)
        b.append(g.bic(x))
    return float(b[0] - b[1])


def main():
    rng = np.random.default_rng(SEED)
    centres = sorted({p.name.split("_")[0] for p in RED.glob("*.parquet")})

    out = {"window_post_onset_days": list(WINDOW), "per_case": [],
           "excluded": {}, }

    cases = []
    for c in centres:
        S = []
        for init in ONSET:
            n, ct = load(c, "nudged", init), load(c, "control", init)
            if n is None or ct is None:
                continue
            S.append(abs(window_means(n, init).mean()
                         - window_means(ct, init).mean()))
        if S and max(S) < 1.0:
            out["excluded"][c] = (f"nudged-minus-control ensemble effect "
                                  f"identically zero (max |S| = {max(S):.3g} Pa)"
                                  f" -- corrupt submission")
            print(f"  EXCLUDED {c}: S = 0 at every init")
            continue
        for init in sorted(ONSET):
            n, ct = load(c, "nudged", init), load(c, "control", init)
            if n is None or ct is None:
                continue
            nm, cm = window_means(n, init), window_means(ct, init)
            if len(nm) < 20 or len(cm) < 20:
                continue
            base, sd = float(cm.mean()), float(cm.std(ddof=1))
            if not np.isfinite(sd) or sd <= 0:
                continue
            # NAM convention: high polar-cap pressure = negative annular mode
            N = -(nm - base) / sd
            C = -(cm - base) / sd
            cases.append({"centre": c, "init": init, "N": N, "C": C,
                          "sd_control_Pa": sd})

    print(f"\n{'centre':8s} {'init':11s} {'nN':>4s} {'nC':>4s} {'sd_ctrl':>8s} "
          f"{'shift':>7s} {'var ratio':>10s} {'KS p':>7s}")
    print("-" * 68)
    for k in cases:
        N, C = k["N"], k["C"]
        shift = float(N.mean() - C.mean())
        vr = float(N.var(ddof=1) / C.var(ddof=1))
        ks = float(stats.ks_2samp(N - N.mean(), C - C.mean()).pvalue)
        k.update(shift=shift, var_ratio=vr, ks_p=ks)
        out["per_case"].append({
            "centre": k["centre"], "init": k["init"],
            "n_nudged": int(len(N)), "n_control": int(len(C)),
            "sd_control_Pa": round(k["sd_control_Pa"], 1),
            "shift_sigma": round(shift, 4),
            "variance_ratio": round(vr, 4),
            "pure_translation_KS_p": round(ks, 4)})
        flag = "  <- short lead" if k["init"] == SHORT_LEAD_INIT else ""
        print(f"{k['centre']:8s} {k['init']:11s} {len(N):4d} {len(C):4d} "
              f"{k['sd_control_Pa']:8.1f} {shift:+7.3f} {vr:10.3f} "
              f"{ks:7.3f}{flag}")

    def pooled(sel, label):
        N = np.concatenate([k["N"] for k in sel])
        C = np.concatenate([k["C"] for k in sel])
        shift = float(N.mean() - C.mean())
        vr = float(N.var(ddof=1) / C.var(ddof=1))
        # bootstrap the variance ratio, resampling MEMBERS within case
        vrs = []
        for _ in range(N_BOOT):
            a = np.concatenate([rng.choice(k["N"], len(k["N"]), replace=True)
                                for k in sel])
            b = np.concatenate([rng.choice(k["C"], len(k["C"]), replace=True)
                                for k in sel])
            vrs.append(a.var(ddof=1) / b.var(ddof=1))
        lo, hi = np.percentile(vrs, [2.5, 97.5])
        ks = stats.ks_2samp(N - N.mean(), C - C.mean())
        # mixture evidence in the forced sample, calibrated on the unforced one
        obs_bic = bic_k2_minus_k1(N - N.mean(), rng)
        null = [bic_k2_minus_k1(rng.choice(C - C.mean(), len(N), replace=True),
                                rng) for _ in range(200)]
        null = np.array([v for v in null if np.isfinite(v)])
        p_mix = float((null >= obs_bic).mean()) if len(null) else np.nan

        # power at a PROPER 5% level. A first version compared the statistic
        # against ONE null draw, which is a coin flip under the null and gave
        # power DECREASING with separation -- impossible, and the giveaway.
        # The critical value is the 95th percentile of the null distribution.
        crit = float(np.percentile(null, 95)) if len(null) else np.nan
        pw = {"critical_value_BIC_at_5pct": round(crit, 3)}
        Cc = C - C.mean()
        for sep in (0.25, 0.5, 0.75, 1.0, 1.5):
            hits = 0
            for _ in range(N_POWER):
                m = rng.random(len(N)) < (2 / 3)
                # two components at the stated separation, mixing fraction 2/3,
                # TOTAL variance held at the observed value so a mixture cannot
                # be detected by spread alone
                x = np.where(m, -sep * (1 - 2 / 3), sep * (2 / 3))
                x = x + rng.choice(Cc, len(N), replace=True)
                x = (x - x.mean()) / x.std(ddof=1) * (N - N.mean()).std(ddof=1)
                hits += (bic_k2_minus_k1(x, rng) > crit)
            pw[f"sep={sep}"] = round(hits / N_POWER, 3)

        res = {"n_cases": len(sel), "n_nudged": int(len(N)),
               "n_control": int(len(C)),
               "shift_sigma": round(shift, 4),
               "variance_ratio": round(vr, 4),
               "variance_ratio_CI95": [round(float(lo), 4), round(float(hi), 4)],
               "forced_variance_sigma2": round(vr - 1.0, 4),
               "pure_translation_KS_stat": round(float(ks.statistic), 4),
               "pure_translation_KS_p": round(float(ks.pvalue), 4),
               "mixture_BIC_p_vs_control": round(p_mix, 4),
               "power_to_detect_mixture": pw}
        out[label] = res
        print(f"\n=== POOLED ({label}) — {len(sel)} cases, "
              f"{len(N)} nudged vs {len(C)} control members ===")
        print(f"  causal shift                : {shift:+.3f} sigma")
        print(f"  variance ratio nudged/control: {vr:.3f} "
              f"[{lo:.3f}, {hi:.3f}]   (1.0 = shape unchanged)")
        print(f"  pure translation, KS        : D = {ks.statistic:.3f}, "
              f"p = {ks.pvalue:.3f}")
        print(f"  two-component mixture       : p = {p_mix:.3f} "
              f"(calibrated on the unforced ensemble)")
        print(f"  power to detect a 2/3 mixture: "
              + ", ".join(f"{k} {v}" for k, v in pw.items()))
        return res

    pooled(cases, "pooled_all")
    keep = [k for k in cases if k["init"] != SHORT_LEAD_INIT]
    if len(keep) < len(cases):
        pooled(keep, "pooled_excluding_short_lead_init")

    (RESULTS / "snapsi_distribution_test.json").write_text(
        json.dumps(out, indent=2), encoding="utf8")
    print("\nSaved -> results/current/8_experiment/snapsi_distribution_test.json")


if __name__ == "__main__":
    main()
