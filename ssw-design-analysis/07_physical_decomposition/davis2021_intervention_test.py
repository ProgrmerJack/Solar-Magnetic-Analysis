#!/usr/bin/env python3
"""
davis2021_intervention_test.py
==============================
A THIRD INDEPENDENT INTERVENTION: DOES REMOVING THE JANUARY 2021 SSW FROM THE
INITIAL STRATOSPHERE CHANGE THE DW-NDW CONTRAST, OR ONLY THE DW RATE?

Registered 2026-10-02, committed BEFORE any member of these ensembles was
downloaded (03_data_ingestion/acquire_davis2022.py writes the reduced data).

DATA (Davis et al. 2022, Nat. Commun. 13, 1136; Zenodo 10.5281/zenodo.5639805,
CC BY 4.0): CESM2(WACCM6) 21-member forecasts, daily, 46 days.
  PRIMARY pair, initialised 4 Jan 2021 (one day before the 5 Jan 2021 SSW):
    SSW     "ssfcst_04jan"   standard initial state           (21/21 produce the SSW)
    no-SSW  "9to12kmRamp"    stratosphere scrambled from an SD run branched
                              16 Nov 2020, troposphere identical  (0/21 produce it)
  SECONDARY pairs, initialised after onset (1 and 8 Feb 2021): "ssfcst_01feb" vs
    "9to12kmRamp_updated" (01feb), "ssfcst_08feb" vs "9to12kmRamp_updated" (08feb);
    there the window is lead days 8-25 after initialisation (the imposed state is
    the post-SSW stratosphere, as for SNAPSI s20190108).
  A different model, a different event and a different kind of intervention
  (initial-state replacement, not nudging) from SNAPSI.

OUTCOME, as close to result L (snapsi_selection_test) as the archive allows
  N = -(polar-cap mean of 1000 hPa height, 60-90N, cos-latitude weights, all
  longitudes) - its mean, divided by its s.d., both taken over the window means of
  the pair's no-SSW arm (the archive has no sea-level pressure; Z1000 is its proxy).
  Window: post-onset days +8..+25 (5 Jan 2021 onset; lead from each file's own
  time axis), daily values. DW = window mean < 0 AND > 50% of days < 0. Contrast
  = mean N(DW) - mean N(NDW), needing >= 3 members per class.
  Regional: northern-Eurasian (50-65N, 10-130E) 2 m temperature, days 8-24, mean
  over the window, standardised by the no-SSW arm.

TESTS (10,000 resamples of members within each arm; seed from the script name)
  R1 rate: DW share, SSW vs no-SSW arm, and the imposed shift s of the mean N.
  R2 PRIMARY: back-shift test as G1 of snapsi_contrast_generality: SSW-arm daily N
     minus s, classified identically; Delta_b = its contrast - no-SSW contrast.
     Also, when both arms form a contrast, the paired difference C_SSW - C_noSSW.
  R3 regional: N-Eurasian T contrast among SSW-arm members against the matched
     shifted null (no-SSW members' N displaced by s, classified, own T kept).
  R4 pooled: R2 and R3 averaged over the three pairs (secondary).
READING, fixed now. With 21 members per arm the intervals will be wide; the
  test can REFUTE contrast invariance but is unlikely to confirm it:
  interval excludes 0 AND |mean| > 0.25 sigma -> "the intervention changed the
  class contrast in this event"; interval inside +-0.25 -> "unchanged within
  tolerance"; otherwise "not rejected; low power". The expected interval
  half-width for n = 21 under one Gaussian population is reported (simulation,
  2,000 draws) so the reader can judge what could have been seen.
PRE-RUN CHECK on synthetic one-population arms (21 x 18 days, shift -0.7 sigma, 200
  draws, before any data): mean Delta_b -0.075, s.d. 0.43; mean resid_T -0.045. So
  this test can only refute changes larger than about 0.8 sigma, and carries a
  small negative finite-sample bias.
CAVEATS stated now: one model; Cohen et al. (2023, Nat. Commun. 14, 3289) argue
  the 9-12 km taper leaves lower-stratospheric influence in the no-SSW arm; the
  polar cap uses 1000 hPa height, not sea-level pressure; 46-day runs.

Output: results/current/8_experiment/davis2021_intervention_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
DATA = HERE.parents[0] / "03_data_ingestion" / "davis2022_reduced.parquet"
NAME = "davis2021_intervention_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, N_SIM, MIN_CLASS, TOL = 10000, 2000, 3, 0.25
ONSET = pd.Timestamp("2021-01-05")
PAIRS = {"04jan": ("ssfcst_04jan", "9to12kmRamp_04jan", "onset"),
         "01feb": ("ssfcst_01feb", "9to12kmRamp_updated_01feb", "init"),
         "08feb": ("ssfcst_08feb", "9to12kmRamp_updated_08feb", "init")}
WIN, TWIN = (8, 25), (8, 24)


def classify(M):
    """M: members x days. Returns window means and DW flags."""
    m = M.mean(1)
    return m, (m < 0) & ((M < 0).mean(1) > 0.5)


def contrast(M, v=None):
    m, dw = classify(M)
    v = m if v is None else v
    if dw.sum() < MIN_CLASS or (~dw).sum() < MIN_CLASS:
        return np.nan
    return float(v[dw].mean() - v[~dw].mean())


def arm_arrays(d, arm, ref):
    """Daily N (members x window days) and window-mean T for one arm."""
    a = d[d.arm == arm]
    t0 = ONSET if ref == "onset" else pd.Timestamp(a.init.iloc[0])
    a = a.assign(day=(pd.to_datetime(a.time) - t0).dt.days)
    w = a[(a.day >= WIN[0]) & (a.day <= WIN[1])].pivot(index="member", columns="day", values="zcap60")
    tw = a[(a.day >= TWIN[0]) & (a.day <= TWIN[1])].groupby("member")["NEURASIA"].mean()
    assert w.shape[1] == WIN[1] - WIN[0] + 1 and not w.isna().any().any(), (arm, w.shape)
    return w, tw


def pair_stats(nud_Z, con_Z, nud_T, con_T):
    """All quantities for one pair; Z in metres, standardised by the no-SSW arm."""
    cm = con_Z.mean(1)
    base, sd = float(cm.mean()), float(cm.std(ddof=1))
    N_s = -(nud_Z.values - base) / sd
    N_c = -(con_Z.values - base) / sd
    mu, sdT = float(con_T.mean()), float(con_T.std(ddof=1))
    T_s = (nud_T.reindex(nud_Z.index).values - mu) / sdT
    T_c = (con_T.reindex(con_Z.index).values - mu) / sdT
    ms, dws = classify(N_s); mc, dwc = classify(N_c)
    s = float(ms.mean() - mc.mean())
    out = {"s": s, "dw_rate_ssw": float(dws.mean()), "dw_rate_nossw": float(dwc.mean()),
           "C_ssw": contrast(N_s), "C_nossw": contrast(N_c), "C_backshift": contrast(N_s - s),
           "CT_ssw": contrast(N_s, T_s), "CT_null": contrast(N_c + s, T_c), "sd_ratio": float(ms.std(ddof=1) / mc.std(ddof=1))}
    out["Delta_b"] = out["C_backshift"] - out["C_nossw"]
    out["paired_diff"] = out["C_ssw"] - out["C_nossw"]
    out["resid_T"] = out["CT_ssw"] - out["CT_null"]
    return out, (N_s, N_c, T_s, T_c)


def boot(N_s, N_c, T_s, T_c, rng):
    keys = ("Delta_b", "paired_diff", "resid_T", "s", "dw_rate_ssw", "dw_rate_nossw")
    bs = {k: [] for k in keys}
    for _ in range(N_BOOT):
        i = rng.integers(len(N_s), size=len(N_s)); j = rng.integers(len(N_c), size=len(N_c))
        Ns, Nc = N_s[i], N_c[j]
        ms, dws = classify(Ns); mc, dwc = classify(Nc)
        s = ms.mean() - mc.mean()
        cb, cc, cs = contrast(Ns - s), contrast(Nc), contrast(Ns)
        bs["Delta_b"].append(cb - cc); bs["paired_diff"].append(cs - cc)
        bs["resid_T"].append(contrast(Ns, T_s[i]) - contrast(Nc + s, T_c[j]))
        bs["s"].append(s); bs["dw_rate_ssw"].append(dws.mean()); bs["dw_rate_nossw"].append(dwc.mean())
    return {k: np.array(v, float) for k, v in bs.items()}


def ci(x):
    x = x[np.isfinite(x)]
    return [round(float(q), 4) for q in np.quantile(x, [0.025, 0.975])] if len(x) > 100 else None


def reading(c, m):
    if c is None or not np.isfinite(m):
        return "not estimable"
    if -TOL <= c[0] and c[1] <= TOL:
        return "unchanged within tolerance"
    if (c[0] > 0 or c[1] < 0) and abs(m) > TOL:
        return "the intervention changed the class contrast in this event"
    return "not rejected; low power"


def expected_halfwidth(rng, n=21, days=18):
    """Interval half-width of Delta_b for two Gaussian arms of n members (one population)."""
    hw = []
    for _ in range(N_SIM // 10):
        def arm(shift):
            return shift + rng.normal(0, 1, (n, 1)) + rng.normal(0, 1.5, (n, days))
        Ns, Nc = arm(-0.5), arm(0.0)
        b = []
        for _ in range(200):
            i = rng.integers(n, size=n); j = rng.integers(n, size=n)
            s = classify(Ns[i])[0].mean() - classify(Nc[j])[0].mean()
            b.append(contrast(Ns[i] - s) - contrast(Nc[j]))
        b = np.array(b)[np.isfinite(b)]
        if len(b) > 50:
            hw.append((np.quantile(b, 0.975) - np.quantile(b, 0.025)) / 2)
    return round(float(np.median(hw)), 3)


def main():
    if not DATA.exists():
        sys.exit("no reduced data: run 03_data_ingestion/acquire_davis2022.py")
    rng = np.random.default_rng(SEED)
    d = pd.read_parquet(DATA)
    res = {"registered": "2026-10-02, committed before any member was downloaded", "seed": SEED,
           "n_boot": N_BOOT, "tolerance_sigma": TOL, "window_days": list(WIN), "t_window_days": list(TWIN),
           "source": "Davis et al. 2022, Zenodo 10.5281/zenodo.5639805 (CC BY 4.0)", "pairs": {}}
    pooled = {k: [] for k in ("Delta_b", "resid_T")}
    for name, (arm_s, arm_c, ref) in PAIRS.items():
        if not {arm_s, arm_c} <= set(d.arm.unique()):
            res["pairs"][name] = "not available"
            continue
        Zs, Ts = arm_arrays(d, arm_s, ref); Zc, Tc = arm_arrays(d, arm_c, ref)
        st, arrs = pair_stats(Zs, Zc, Ts, Tc)
        b = boot(*arrs, rng)
        blk = {k: round(v, 4) if isinstance(v, float) else v for k, v in st.items()}
        blk["n_members"] = [int(len(Zs)), int(len(Zc))]
        blk["ci95"] = {k: ci(v) for k, v in b.items()}
        blk["reading_Delta_b"] = reading(blk["ci95"]["Delta_b"], st["Delta_b"])
        blk["reading_resid_T"] = reading(blk["ci95"]["resid_T"], st["resid_T"])
        res["pairs"][name] = blk
        for k in pooled:
            pooled[k].append(b[k])
        print(name, {k: blk[k] for k in ("s", "dw_rate_ssw", "dw_rate_nossw", "C_ssw", "C_nossw", "Delta_b",
                                          "paired_diff", "resid_T", "sd_ratio")}, blk["ci95"]["Delta_b"],
              blk["reading_Delta_b"], "| T", blk["ci95"]["resid_T"], blk["reading_resid_T"], flush=True)
    if len(pooled["Delta_b"]) > 1:
        res["R4_pooled_mean_of_pairs"] = {}
        for k, v in pooled.items():
            n = min(len(x) for x in v)
            m = np.nanmean(np.vstack([x[:n] for x in v]), axis=0)
            pt = float(np.nanmean([res["pairs"][p][k] for p in res["pairs"] if isinstance(res["pairs"][p], dict)]))
            res["R4_pooled_mean_of_pairs"][k] = {"mean": round(pt, 4), "ci95": ci(m), "reading": reading(ci(m), pt)}
        print("pooled", res["R4_pooled_mean_of_pairs"])
    res["expected_halfwidth_one_population_n21"] = expected_halfwidth(rng)
    print("expected half-width (one population, n = 21):", res["expected_halfwidth_one_population_n21"])
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
