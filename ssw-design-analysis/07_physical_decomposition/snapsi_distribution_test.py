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
  - s20190108 initialises 2019-01-08 with onset 2019-01-02, so onset PRECEDES
    the initialisation by 6 days and the +8..+25 post-onset window falls at lead
    2-19 days, before the ensemble has fully diverged; its control spread is
    smaller than elsewhere. The headline pool (`pooled_all`) INCLUDES it; the
    pool without it is reported alongside as the sensitivity.
  - NRL is excluded: its nudged-minus-control effect is zero to within 0.023 Pa
    per member at three of four initialisations (defective submission).

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
        # member values in control-sd units (NAM sign), for the display item
        out.setdefault("member_A", {})[f"{k['centre']}|{k['init']}"] = {
            "nudged": [round(float(x), 4) for x in N],
            "control": [round(float(x), 4) for x in C]}
        flag = "  <- short lead" if k["init"] == SHORT_LEAD_INIT else ""
        print(f"{k['centre']:8s} {k['init']:11s} {len(N):4d} {len(C):4d} "
              f"{k['sd_control_Pa']:8.1f} {shift:+7.3f} {vr:10.3f} "
              f"{ks:7.3f}{flag}")

    def pooled(sel, label):
        # EACH ENSEMBLE IS RE-CENTRED ON ITS OWN MEAN BEFORE POOLING.
        #
        # Standardisation uses the CONTROL mean and sd of each (centre, init),
        # so the control arm lands at ~0 in every case while the nudged arm
        # lands on THAT case's causal shift. Concatenating without re-centring
        # therefore adds the between-ensemble spread of those shifts to the
        # nudged variance and to nothing else, and the shape tests then measure
        # that spread rather than the within-ensemble dispersion they are for.
        #
        # It stayed invisible at 3 centres because the shift spread was only
        # 0.295 sigma, contributing 0.087 to a ratio of 0.941. At 9 centres the
        # spread is 1.07 sigma, contributing 1.15 -- more than the quantity
        # being estimated -- and the uncorrected ratio reads 2.085 while the
        # within-ensemble ratio is 0.952. The correction does not overturn the
        # published conclusion, it is what lets it survive the larger sample.
        #
        # The shift itself is still measured, between ensembles, where it
        # belongs -- and reported with its spread, which is a real quantity.
        Nc_by_case = [k["N"] - k["N"].mean() for k in sel]
        Cc_by_case = [k["C"] - k["C"].mean() for k in sel]
        N = np.concatenate(Nc_by_case)
        C = np.concatenate(Cc_by_case)
        case_shifts = np.array([float(k["N"].mean() - k["C"].mean()) for k in sel])
        shift = float(case_shifts.mean())
        vr = float(N.var(ddof=1) / C.var(ddof=1))
        # bootstrap the variance ratio, resampling MEMBERS within case
        vrs = []
        for _ in range(N_BOOT):
            a = np.concatenate([rng.choice(x, len(x), replace=True)
                                for x in Nc_by_case])
            b = np.concatenate([rng.choice(x, len(x), replace=True)
                                for x in Cc_by_case])
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
               # Between-ensemble spread of the causal shift. NOT sigma_f: it
               # mixes between-MODEL with between-EVENT variation, since the
               # ensembles span 9 centres and only 2 events. Reported because
               # it is what contaminates an uncentred pooled variance ratio,
               # and because it is large -- it must not be quoted as an
               # event-to-event forced spread.
               "shift_sd_across_cases": round(float(case_shifts.std(ddof=1)), 4),
               "shift_range_across_cases": [round(float(case_shifts.min()), 4),
                                            round(float(case_shifts.max()), 4)],
               "variance_ratio_is_within_ensemble": True,
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
        print(f"  causal shift                : {shift:+.3f} sigma "
              f"(sd across cases {case_shifts.std(ddof=1):.3f}, "
              f"range {case_shifts.min():+.2f}..{case_shifts.max():+.2f})")
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
        json.dumps(out, indent=2), encoding="utf8", newline="\n")
    print("\nSaved -> results/current/8_experiment/snapsi_distribution_test.json")


if __name__ == "__main__":
    main()
