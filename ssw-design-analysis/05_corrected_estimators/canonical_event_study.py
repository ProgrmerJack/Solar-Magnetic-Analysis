#!/usr/bin/env python3
"""
canonical_event_study.py
========================
The canonical AO/NAO lead-lag profile, computed on the FROZEN catalogue and
repeated across every frozen sensitivity set.

This serves two purposes:
  * it is the physical anchor of the paper (Figure 3 in the plan);
  * it is Gate 5 condition 7 -- the profile shape must not be qualitatively
    dependent on which catalogue is used.

WHY IT MUST BE RECOMPUTED
  The superseded event studies (r113/r114) used a 47-event list, which was the
  UNION across reanalyses. The preregistered primary catalogue is the >=4/6
  consensus rule, giving 39 events in 33 winters. Every number carried forward
  from the old runs is therefore keyed to the wrong exposure and must be
  re-derived here.

SPECIFICATION (all fixed before running)
  seasonal model   3 annual harmonics -- the Gate 3 primary (smallest bias)
  bins             mutually exclusive, nearest-onset assignment
  baseline         winter days >75 d from ANY onset in the set being tested
  inference        winter-block bootstrap, 1,200 replicates (Gate 3 floor: 1,000)
  reported         the WHOLE profile, both N_events and N_unique_winters

Nothing is selected on the outcome. All eight catalogues are reported whatever
they show.

Output: canonical_event_study.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "2_event_study"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue                    # noqa: E402
import calibrate_seasonality as C                             # noqa: E402

OUT = RESULTS / "canonical_event_study.json"   # re-run under frozen Gate 3 rules

BINS = [(-60, -46), (-45, -31), (-30, -16), (-15, -1),
        (0, 14), (15, 29), (30, 44), (45, 60)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
BASELINE_GAP = 75
SPEC = "harm3"          # Gate 3 primary: smallest residual bias (+0.033), best calibration
N_BOOT = 1200          # Gate 3 standing rule: never below 1,000
SETS = ["primary", "primary_compendium_only", "consensus_strict", "consensus_half", "union",
        "era5", "jra_55", "ncep_ncar", "merra2"]


def load_series(name):
    f = {"AO": "ao_daily_cpc.txt", "NAO": "nao_daily_cpc.txt"}[name]
    s = C.read_cpc(ROOT / "data/processed/atmospheric" / f)
    d = pd.DataFrame({"y": s.values}, index=s.index)
    d = d[np.isin(d.index.month, C.SEASON)]
    d["winter"] = C.winter_of(d.index)
    d["doy"] = d.index.dayofyear
    return d


def build(d, onsets):
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets]))
    dd = d.index.values.astype("datetime64[D]")
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    code = np.full(len(dd), -1)
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    out = d.copy()
    out["bin"] = code
    out["is_baseline"] = np.abs(lag) > BASELINE_GAP
    return out[(out["bin"] >= 0) | out["is_baseline"]].copy()


def fit(dd):
    y = dd["y"].values.astype(float)
    S = C.harmonics(dd["doy"].values, 3)   # SPEC = harm3
    D = np.column_stack([(dd["bin"].values == i).astype(float)
                         for i in range(len(BINS))])
    A = np.column_stack([y, D, S])
    Ad = C_demean(A, dd["winter"].values)
    try:
        beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return beta[:len(BINS)]
    except Exception:
        return None


def C_demean(A, codes):
    _, c = np.unique(codes, return_inverse=True)
    n = c.max() + 1
    cnt = np.bincount(c, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (s / cnt)[c]
    return out


def run(d, onsets, seed=0):
    dd = build(d, onsets)
    est = fit(dd)
    if est is None:
        return None
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(dd["winter"].unique()))
    idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        rows = np.concatenate([idx[w] for w in pick])
        sub = dd.iloc[rows].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        b = fit(sub)
        if b is not None and np.all(np.isfinite(b)):
            boot.append(b)
    boot = np.array(boot)
    prof = {}
    for i, lab in enumerate(LABELS):
        col = boot[:, i]
        prof[lab] = {
            "effect": round(float(est[i]), 4),
            "CI95": [round(float(np.percentile(col, 2.5)), 4),
                     round(float(np.percentile(col, 97.5)), 4)],
            "p_two_sided": float(2 * min((col >= 0).mean(), (col <= 0).mean())),
        }
    return prof, dd


def main():
    res = {"seasonal_spec": SPEC, "bins": LABELS,
           "baseline_gap_days": BASELINE_GAP, "n_boot": N_BOOT, "outcomes": {}}

    for name in ("AO", "NAO"):
        d = load_series(name)
        res["outcomes"][name] = {}
        print(f"\n{'='*78}\n### {name}  ({d.index.min().date()}..{d.index.max().date()}, "
              f"{d['winter'].nunique()} winters)")
        for s in SETS:
            onsets = load_catalogue(s)
            onsets = onsets[(onsets >= d.index.min()) & (onsets <= d.index.max())]
            if len(onsets) < 10:
                print(f"  {s:26s} skipped ({len(onsets)} events in range)")
                continue
            # crc32, not hash(): Python randomises string hashing per process
            # (PYTHONHASHSEED), so hash(s) drew a different bootstrap sample on
            # every run and no p-value here was reproducible. See seed_stability.py.
            out = run(d, onsets, seed=zlib.crc32(s.encode()) % 10000)
            if out is None:
                continue
            prof, dd = out
            n_ev = len(onsets)
            n_w = int(pd.DatetimeIndex(onsets).to_series()
                      .groupby(C.winter_of(pd.DatetimeIndex(onsets))).ngroups)
            res["outcomes"][name][s] = {
                "n_events": n_ev, "n_unique_winters": n_w,
                "n_days": int(len(dd)),
                "n_baseline_days": int(dd["is_baseline"].sum()),
                "profile": prof,
            }
            star = "".join("*" if prof[l]["p_two_sided"] < 0.05 else "." for l in LABELS)
            vals = " ".join(f"{prof[l]['effect']:+.2f}" for l in LABELS)
            print(f"  {s:26s} n={n_ev:3d}/{n_w:2d}w  {vals}   {star}")

        # Gate 5 condition 7: is the shape stable across catalogues?
        sets_done = list(res["outcomes"][name])
        post = [i for i, l in enumerate(LABELS) if l.startswith("+")]
        pre = [i for i, l in enumerate(LABELS) if l.startswith("-")]
        signs_post, signs_pre = [], []
        for s in sets_done:
            p = res["outcomes"][name][s]["profile"]
            signs_post.append(all(p[LABELS[i]]["effect"] < 0 for i in post))
            signs_pre.append(all(p[LABELS[i]]["effect"] < 0 for i in pre))
        res["outcomes"][name]["_stability"] = {
            "catalogues_tested": len(sets_done),
            "all_post_onset_bins_negative_in_every_catalogue": bool(all(signs_post)),
            "all_pre_onset_bins_negative_in_every_catalogue": bool(all(signs_pre)),
        }
        print(f"  -> post-onset negative in every catalogue: {all(signs_post)}; "
              f"pre-onset negative in every catalogue: {all(signs_pre)}")

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
