#!/usr/bin/env python3
"""
gate3_clean_null.py
===================
Corrects a flaw found in the previous Gate 3 attempt.

WHAT WENT WRONG
  `finalise_gate3.py` implemented the requested "exclude pseudo-onsets from the
  real-event influence zone" check. It produced:

      full_pool       null = -0.006   (essentially zero)
      excluded        null = +0.234   coverage 0.878, FPR 0.122
      nonssw_winters  null = +0.149

  The three draws were required to agree; they do not, and `excluded` is biased
  POSITIVE by about a quarter of a standard deviation.

  The cause is the baseline, not the onsets. In an event study the baseline is
  "days far from any onset" -- but under that draw, "onset" means the FAKE onset.
  Pushing fake onsets away from real events therefore leaves the real,
  genuinely depressed SSW periods sitting inside the BASELINE. The fake bins are
  then compared against a reference contaminated with the real negative
  response, so they appear positive. +0.234 is the mirror image of the real
  effect leaking into the reference category.

  This is the same class of error as the original single-window flaw: the
  comparison group, not the treatment group, was contaminated.

THE FIX
  Remove real-event influence from the DATA before drawing anything, so that
  neither the fake bins nor the fake baseline can contain it:

    1. drop every day within [-60, +75] d of any REAL onset;
    2. draw pseudo-onsets in the surviving clean days, preserving the observed
       day-of-year distribution;
    3. fit the event study on the cleaned series.

  Now every day in the calibration -- treated and control alike -- is free of
  real SSW influence, and a non-zero result can only come from the estimator or
  the seasonal model.

Output: gate3_clean_null.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import calibrate_seasonality as C

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_event_catalogues"))
from build_catalogue import load_catalogue      # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
OUT = RESULTS / "gate3_clean_null.json"
CKPT = RESULTS / "_gate3_clean_null_checkpoint.json"
CKPT_EVERY = 25          # save progress this often during the coverage loop


def load_ckpt():
    return json.loads(CKPT.read_text(encoding="utf8")) if CKPT.exists() else {}


def save_ckpt(ck):
    CKPT.write_text(json.dumps(ck, indent=1), encoding="utf8")

N_MEAN = 4000
N_COVER = 300   # as in GATE3_CLOSURE.md
N_BOOT = 1000   # STANDING RULE (GATE3_CLOSURE.md): no reported interval may
                # use fewer than 1000. This script shipped at 150, which is
                # the setting that rule forbids -- and every artifact it wrote
                # says passes=False, while GATE3_CLOSURE.md records Gate 3 as
                # CLOSED on 0.947/0.053 obtained at 1000. Nothing on disk
                # supported the closure until this was corrected.
INFLUENCE = (-60, 75)        # real-event zone removed from the data entirely
SPECS = ["harm3", "harm6", "cyclic_spline12"]   # PRIMARY first: the
    # closure rests on harm3, so compute it before the sensitivity specs


def clopper_pearson(k, n, alpha=0.05):
    lo = stats.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return round(float(lo), 3), round(float(hi), 3)


def real_influence_mask(index, real_onsets):
    """True where a day lies within the influence zone of ANY real onset."""
    dd = pd.DatetimeIndex(index).values.astype("datetime64[D]")
    m = np.zeros(len(dd), bool)
    for o in real_onsets:
        lag = (dd - np.datetime64(pd.Timestamp(o), "D")).astype(int)
        m |= (lag >= INFLUENCE[0]) & (lag <= INFLUENCE[1])
    return m


def build_by_doy(clean_index):
    """Day-of-year -> candidate days. Built ONCE; see draw_clean."""
    by_doy = {}
    for t in clean_index:
        by_doy.setdefault(t.dayofyear, []).append(t)
    return by_doy


def draw_clean(clean_index, doys, rng, by_doy=None):
    """Pseudo-onsets on clean days only, preserving the observed DOY set.

    `by_doy` is hoisted out because this used to rebuild a dict over all 8,862
    clean timestamps on EVERY call -- 4,000 times in the null-mean loop alone,
    per spec. Passing it in changes no RNG call and so leaves results
    bit-identical; it only stops the same index being rebuilt 12,000 times.
    """
    if by_doy is None:
        by_doy = build_by_doy(clean_index)
    out = []
    for doy in rng.permutation(doys):
        cand = by_doy.get(int(doy))
        if cand:
            out.append(cand[rng.integers(len(cand))])
    return pd.DatetimeIndex(out)


def main():
    ao = C.read_cpc(C.ROOT / "data/processed/atmospheric/ao_daily_cpc.txt")
    d = pd.DataFrame({"y": ao.values}, index=ao.index)
    d = d[np.isin(d.index.month, C.SEASON)]
    d["winter"] = C.winter_of(d.index)
    d["doy"] = d.index.dayofyear

    onsets = load_catalogue("primary")
    onsets = onsets[(onsets >= d.index.min()) & (onsets <= d.index.max())]

    infl = real_influence_mask(d.index, onsets)
    clean = d[~infl].copy()
    doys = pd.DatetimeIndex(onsets).dayofyear.values
    print(f"AO winter days {len(d):,}; removed {infl.sum():,} within "
          f"[{INFLUENCE[0]:+d},{INFLUENCE[1]:+d}] d of a real onset "
          f"({100*infl.mean():.0f}%); {len(clean):,} clean days remain "
          f"across {clean['winter'].nunique()} winters")

    # confirm the diagnosis: how contaminated was the old 'excluded' baseline?
    diag = {}
    rng = np.random.default_rng(5)
    fake = draw_clean(clean.index, doys, rng)
    codes, lag = None, None
    from calibrate_seasonality import assign_bins
    code, lg = assign_bins(d.index, fake)
    base = (code < 0) & (np.abs(lg) > C.BASELINE_GAP)
    diag["old_design_baseline_contamination"] = round(
        float(infl[base].mean()), 3)
    print(f"diagnosis: under the old 'excluded' draw, "
          f"{100*diag['old_design_baseline_contamination']:.0f}% of baseline days "
          f"lay inside real-event influence -> the positive bias")

    res = {"n_mean_draws": N_MEAN, "n_cover_draws": N_COVER,
           "influence_zone": list(INFLUENCE),
           "n_days_total": int(len(d)), "n_days_clean": int(len(clean)),
           "diagnosis": diag, "specs": {}}

    by_doy = build_by_doy(clean.index)          # built once, reused everywhere

    for spec in SPECS:
        # the null-mean loop is deterministic (seed 11) and costs 4,000 fits,
        # so it is cached too -- otherwise every resume recomputes it for every
        # already-finished spec before reaching the work that is left
        ck = load_ckpt()
        mkey = f"{spec}|MEAN|N_MEAN={N_MEAN}"
        if mkey in ck:
            ests = np.array(ck[mkey])
            print(f"  {spec:16s} null-mean loop from cache ({len(ests)} draws)")
        else:
            rng = np.random.default_rng(11)
            ests = []
            for _ in range(N_MEAN):
                fk = draw_clean(clean.index, doys, rng, by_doy)
                if len(fk) < 5:
                    continue
                b = C.fit_once(clean, fk, spec)
                if b is not None and np.isfinite(b[C.TEST_BIN]):
                    ests.append(b[C.TEST_BIN])
            ck[mkey] = [float(x) for x in ests]
            save_ckpt(ck)
            ests = np.array(ests)

        # ---- coverage, CHECKPOINTED --------------------------------------
        # At N_BOOT=1000 this is ~20 min per spec, longer than a single run
        # slot. Progress is saved every CKPT_EVERY draws and resumed on the
        # next invocation, so the job completes across several runs and a kill
        # never costs more than CKPT_EVERY draws.
        ck = load_ckpt()
        key = f"{spec}|N_COVER={N_COVER}|N_BOOT={N_BOOT}"
        st = ck.get(key, {"done": 0, "cov": 0, "rej": 0})
        rng = np.random.default_rng(23)
        # replay the generator so a resumed run continues the same draw sequence
        for _ in range(st["done"]):
            draw_clean(clean.index, doys, rng, by_doy)
            rng.integers(1e6)
        cov, rej = st["cov"], st["rej"]
        if st["done"]:
            print(f"  {spec:16s} resuming coverage at draw {st['done']}/{N_COVER}")
        for i in range(st["done"], N_COVER):
            fk = draw_clean(clean.index, doys, rng, by_doy)
            seed = int(rng.integers(1e6))
            if len(fk) >= 5:
                ci = C.bootstrap_ci(clean, fk, spec, n_boot=N_BOOT, seed=seed)
                if ci:
                    if ci[0] <= 0 <= ci[1]:
                        cov += 1
                    else:
                        rej += 1
            if (i + 1) % CKPT_EVERY == 0 or i + 1 == N_COVER:
                ck[key] = {"done": i + 1, "cov": cov, "rej": rej}
                save_ckpt(ck)
                print(f"    {spec}: {i+1}/{N_COVER} draws, "
                      f"cov={cov/(cov+rej):.3f}" if cov + rej else "")
        n = cov + rej
        e = {
            "null_mean": round(float(ests.mean()), 4),
            "null_sd": round(float(ests.std()), 4),
            "mean_CI95": [round(float(ests.mean() - 1.96 * ests.std() / np.sqrt(len(ests))), 4),
                          round(float(ests.mean() + 1.96 * ests.std() / np.sqrt(len(ests))), 4)],
            "n_draws": int(len(ests)),
            "coverage": round(cov / n, 3) if n else None,
            "coverage_CI95": clopper_pearson(cov, n) if n else None,
            "fpr": round(rej / n, 3) if n else None,
            "fpr_CI95": clopper_pearson(rej, n) if n else None,
            "n_cover": int(n),
        }
        e["passes"] = bool(
            abs(e["null_mean"]) <= 0.05
            and e["coverage_CI95"] and e["coverage_CI95"][0] <= 0.95 <= e["coverage_CI95"][1]
            and e["fpr_CI95"] and e["fpr_CI95"][0] <= 0.05 <= e["fpr_CI95"][1])
        res["specs"][spec] = e
        res["n_boot"] = N_BOOT
        res["partial"] = len(res["specs"]) < len(SPECS)
        OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
        print(f"  {spec:16s} null={e['null_mean']:+.4f} {e['mean_CI95']}  "
              f"cov={e['coverage']} {e['coverage_CI95']}  "
              f"FPR={e['fpr']} {e['fpr_CI95']}  pass={e['passes']}")

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\npassing:", [s for s, v in res["specs"].items() if v["passes"]] or "NONE")
    print(f"Saved -> {OUT}")


if __name__ == "__main__":
    main()
