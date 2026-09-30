#!/usr/bin/env python3
"""
finalise_gate3.py
=================
FLAWED. DO NOT CITE. SUPERSEDED BY `gate3_clean_null.py`.

Its `excluded` pseudo-onset draw leaves real-event-influenced days inside the
baseline, so the null centres near +0.29 sigma instead of zero
(`FINDING_gate3_contamination.md` line 11, which labels this script "flawed").
The name is the trap: it implies this is the final word on Gate 3 and it is not.
Gate 3 closes on `gate3_clean_null.py` -- 3 harmonics, null mean +0.0052
[-0.0062, +0.0166], coverage 0.930, FPR 0.070, at N_BOOT=1000. Its output
`gate3_final.json` was moved to `results/superseded/flawed_producer/` on
2026-09-17; running this script would write it back.

Kept, not deleted, because the contamination it exposed is the evidence for the
corrected design. Read it; do not run it for a result.

--- original header follows ---

Closes Gate 3 with the two checks required by review.

CHECK 1 — pseudo-onsets must not sit inside real-event influence.
  Drawing from all 77 winters removed the gross contamination, but a fake onset
  can still land near a real SSW inside an SSW-hosting winter. Three draws are
  compared; they must agree within Monte Carlo error:

    full_pool       any winter, observed day-of-year distribution (previous test)
    excluded        as above, but rejecting any date within [-60, +75] d of ANY
                    real onset -- no overlap with real lead-lag influence
    nonssw_winters  drawn only from winters containing no SSW at all

CHECK 2 — enough draws to freeze a false-positive rate.
  400 draws puts ~20 events behind an FPR of 0.05. Null means use N_MEAN draws;
  coverage/FPR use N_COVER draws each with its own winter-block bootstrap
  interval, and both are reported with exact binomial (Clopper-Pearson) intervals.

Specifications carried forward: the two clean passes (6 harmonics primary,
12-knot cyclic spline sensitivity) plus 3 harmonics for continuity with the
superseded work.

Output: gate3_final.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import calibrate_seasonality as C

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_event_catalogues"))
from build_catalogue import load_catalogue      # the ONLY permitted event source

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
OUT = RESULTS / "gate3_final.json"

N_MEAN = 4000        # draws for the null mean (as in GATE3_CLOSURE.md)
N_COVER = 300        # draws that also get a bootstrap interval
COVER_KINDS = ("excluded",)   # coverage on the strictest draw only
N_BOOT = 1000        # STANDING RULE (GATE3_CLOSURE.md): no reported interval
                     # may use <1000. At 150 the same estimator gives coverage
                     # 0.920 / FPR 0.080 instead of 0.947 / 0.053, so the 150
                     # this script previously used contradicted its own closure.
EXCLUDE = (-60, 75)  # real-event influence zone to keep pseudo-onsets out of
SPECS = ["harm6", "cyclic_spline12", "harm3"]


def clopper_pearson(k, n, alpha=0.05):
    lo = stats.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return round(float(lo), 3), round(float(hi), 3)


def draw_pseudo(kind, real_onsets, all_w, nonssw_w, real_days, rng):
    """Pseudo-onsets preserving the observed day-of-year distribution."""
    doys = pd.DatetimeIndex(real_onsets).dayofyear.values
    pool = nonssw_w if kind == "nonssw_winters" else all_w
    if len(pool) == 0:
        return pd.DatetimeIndex([])
    out = []
    for doy in rng.permutation(doys):
        for _ in range(40):                       # rejection sampling
            w = int(rng.choice(pool))
            year = w - 1 if doy > 250 else w
            try:
                t = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=int(doy) - 1)
            except Exception:
                continue
            if kind == "excluded":
                d = (real_days - np.datetime64(t, "D")).astype(int)
                # reject if this date lies in any real event's influence zone
                if np.any((-d >= EXCLUDE[0]) & (-d <= EXCLUDE[1])):
                    continue
            out.append(t)
            break
    return pd.DatetimeIndex(out)


def main():
    ao = C.read_cpc(C.ROOT / "data/processed/atmospheric/ao_daily_cpc.txt")
    d = pd.DataFrame({"y": ao.values}, index=ao.index)
    d = d[np.isin(d.index.month, C.SEASON)]
    d["winter"] = C.winter_of(d.index)
    d["doy"] = d.index.dayofyear

    # frozen primary catalogue: consensus >=2/3 of COVERING products,
    # 43 events / 36 winters after the 2026-07-30 re-freeze
    onsets = load_catalogue("primary")
    onsets = onsets[(onsets >= d.index.min()) & (onsets <= d.index.max())]
    real_days = np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets])

    all_w = np.array(sorted(d["winter"].unique()))
    ssw_w = set(C.winter_of(pd.DatetimeIndex(onsets)))
    nonssw_w = np.array([w for w in all_w if w not in ssw_w])
    print(f"AO: {len(d):,} winter days, {len(all_w)} winters "
          f"({len(nonssw_w)} with no SSW), {len(onsets)} events")

    res = {"n_mean_draws": N_MEAN, "n_cover_draws": N_COVER,
           "exclusion_zone_days": list(EXCLUDE),
           "n_winters_total": int(len(all_w)),
           "n_winters_nonssw": int(len(nonssw_w)), "specs": {}}

    for spec in SPECS:
        res["specs"][spec] = {}
        for kind in ("full_pool", "excluded", "nonssw_winters"):
            rng = np.random.default_rng(11)
            ests = []
            for _ in range(N_MEAN):
                fake = draw_pseudo(kind, onsets, all_w, nonssw_w, real_days, rng)
                if len(fake) < 5:
                    continue
                b = C.fit_once(d, fake, spec)
                if b is not None and np.isfinite(b[C.TEST_BIN]):
                    ests.append(b[C.TEST_BIN])
            ests = np.array(ests)

            rng = np.random.default_rng(23)
            cov = rej = 0
            n_cover_draws = N_COVER if kind in COVER_KINDS else 0
            for _ in range(n_cover_draws):
                fake = draw_pseudo(kind, onsets, all_w, nonssw_w, real_days, rng)
                if len(fake) < 5:
                    continue
                ci = C.bootstrap_ci(d, fake, spec, n_boot=N_BOOT, seed=int(rng.integers(1e6)))
                if ci:
                    if ci[0] <= 0 <= ci[1]:
                        cov += 1
                    else:
                        rej += 1
            n = cov + rej
            entry = {
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
            res["specs"][spec][kind] = entry
            print(f"  {spec:16s} {kind:15s} null={entry['null_mean']:+.4f} "
                  f"{entry['mean_CI95']}  cov={entry['coverage']} {entry['coverage_CI95']}  "
                  f"FPR={entry['fpr']} {entry['fpr_CI95']}")

    # agreement between draws, and pass decision on the strictest (excluded)
    for spec, v in res["specs"].items():
        means = [v[k]["null_mean"] for k in v]
        v_ = v["excluded"]
        res["specs"][spec]["agreement_range"] = round(max(means) - min(means), 4)
        res["specs"][spec]["passes_strict"] = bool(
            abs(v_["null_mean"]) <= 0.05
            and v_["coverage_CI95"][0] <= 0.95 <= v_["coverage_CI95"][1]
            and v_["fpr_CI95"][0] <= 0.05 <= v_["fpr_CI95"][1])

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\npasses under the contamination-excluded draw:",
          [s for s, v in res["specs"].items() if v["passes_strict"]] or "NONE")
    print(f"Saved -> {OUT}")


if __name__ == "__main__":
    main()
