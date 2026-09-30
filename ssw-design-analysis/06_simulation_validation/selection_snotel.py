#!/usr/bin/env python3
"""
selection_snotel.py
===================
Extends selection_on_outcome.py from circulation indices to a REAL surface
impact variable.

WHY THIS MATTERS MORE THAN THE INDEX VERSION
  Showing that selecting on the AO biases the AO is nearly tautological. The
  claim that bites is that the bias propagates into the applied impact variables
  studies actually care about -- snowpack, temperature, air quality, energy
  demand -- because those are driven by the same surface circulation the
  selection conditioned on.

  Here the outcome is the SNOTEL network: 945 stations, 10.6M daily
  observations, western US. Mean daily air temperature and snow depth.

  Note the geography honestly: SNOTEL is western US, where the AO/NAO signal is
  weak and the Pacific (PNA) pathway dominates. So the SSW effect itself may be
  small or absent here. That does NOT weaken the test -- it sharpens it. If a
  classification built from the ATLANTIC-sector surface NAM manufactures an
  apparent western-US snowpack impact out of a true null, that is the cleanest
  possible demonstration of selection bias, because there is little real signal
  for it to hide behind.

DESIGN
  Identical to selection_on_outcome.py:
    - events classified dSSW by the published Karpechko surface-NAM criteria
    - pseudo-onsets drawn on CLEAN days (real-event zone deleted first)
    - day-of-year climatology built from CLEAN days only
    - null selection RATE-MATCHED to the real classification rate
  so any apparent impact in the null is pure selection.

Output: selection_snotel.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "4_selection_bias"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as S                    # noqa: E402

ROOT = HERE.parents[1]
SNOTEL = ROOT / "data" / "processed" / "cryosphere" / "snotel_daily.parquet"
N_DRAW = 800
VARS = {"tavg_c": "mean daily air temperature", "snwd_mm": "snow depth"}


def load_snotel():
    """Network daily mean per variable, standardised, winter days only."""
    d = pd.read_parquet(SNOTEL, columns=["station_id", "tavg_c", "snwd_mm"])
    d.index = pd.to_datetime(d.index)
    # snotel_daily.parquet is tz-aware (UTC) while the catalogue and the CPC
    # indices are tz-naive; comparing the two raises rather than silently
    # misaligning, so strip the zone here at the single point of entry
    if d.index.tz is not None:
        d.index = d.index.tz_convert("UTC").tz_localize(None)
    out = {}
    for v in VARS:
        s = d.groupby(d.index)[v].mean().sort_index()
        s = s[np.isin(s.index.month, M.SEASON)] if hasattr(M, "SEASON") else s
        s = s[s.index.month.isin((11, 12, 1, 2, 3, 4))]
        s = (s - s.mean()) / s.std()          # standardise so units match indices
        out[v] = s.dropna()
    return out


def main():
    ao = M.load("ao")["y"]
    snam = S.strat_nam()
    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]

    sno = load_snotel()
    for v, s in sno.items():
        print(f"SNOTEL {v}: {len(s):,} winter days "
              f"{s.index.min().date()}..{s.index.max().date()}")
    # restrict events to the SNOTEL period
    lo = max(s.index.min() for s in sno.values())
    hi = min(s.index.max() for s in sno.values())
    real_s = real[(real >= lo) & (real <= hi)]
    print(f"events inside SNOTEL period: {len(real_s)} of {len(real)}")

    lab = S.classify(real_s, ao, snam)
    rate_real = float(lab.mean())
    print(f"  dSSW {int(lab.sum())}/{len(real_s)} ({100 * rate_real:.0f}%)")

    mask_full = G.real_influence_mask(ao.index, real)
    clean_idx = G.zone_free_index(ao.index, real)
    clean_idx = clean_idx[(clean_idx >= lo) & (clean_idx <= hi)]
    clims = {}
    for v, s in sno.items():
        cs = s[~G.real_influence_mask(s.index, real)]
        clims[v] = cs.groupby(cs.index.dayofyear).mean()

    real_all, real_d = {}, {}
    print(f"\n{'variable':10s} {'ALL events':>12s} {'dSSW only':>12s}")
    print("-" * 36)
    for v, s in sno.items():
        pe = S.per_event(s, real_s, clims[v])
        real_all[v] = float(np.nanmean(pe))
        real_d[v] = float(np.nanmean(pe[lab]))
        print(f"{v:10s} {real_all[v]:+12.3f} {real_d[v]:+12.3f}")

    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real_s)])
    rng = np.random.default_rng(20260731)
    null = {v: {"all": [], "matched": []} for v in sno}
    for i in range(N_DRAW):
        fake = G.draw_clean(clean_idx, doys, rng)
        if len(fake) < 8:
            continue
        sel_ao = S.per_event(ao, fake, ao.groupby(ao.index.dayofyear).mean())
        k = max(1, int(round(rate_real * len(fake))))
        order = np.argsort(np.where(np.isnan(sel_ao), np.inf, sel_ao))
        matched = np.zeros(len(fake), bool)
        matched[order[:k]] = True
        for v, s in sno.items():
            pe = S.per_event(s, fake, clims[v])
            if np.all(np.isnan(pe)):
                continue
            null[v]["all"].append(np.nanmean(pe))
            null[v]["matched"].append(np.nanmean(pe[matched]))
        if (i + 1) % 200 == 0:
            print(f"  {i + 1}/{N_DRAW} null draws", flush=True)

    res = {"n_draw": N_DRAW, "n_events": int(len(real_s)),
           "rate_real": round(rate_real, 3), "variables": {}}
    print(f"\n{'variable':10s} {'real dSSW':>11s} {'null all':>10s} "
          f"{'null matched':>13s} {'reproduced':>11s}")
    print("-" * 60)
    for v in sno:
        a = np.array(null[v]["all"])
        m = np.array(null[v]["matched"])
        rep = 100 * m.mean() / real_d[v] if real_d[v] else np.nan
        res["variables"][v] = {
            "description": VARS[v],
            "real_all_events": round(real_all[v], 4),
            "real_dSSW_only": round(real_d[v], 4),
            "null_all_events": round(float(a.mean()), 4),
            "null_rate_matched": round(float(m.mean()), 4),
            "null_matched_CI95": [round(float(np.percentile(m, 2.5)), 4),
                                  round(float(np.percentile(m, 97.5)), 4)],
            "pct_reproduced_by_selection": round(float(rep), 1),
            "p_real_vs_matched_null": float((m <= real_d[v]).mean())
            if real_d[v] < 0 else float((m >= real_d[v]).mean())}
        print(f"{v:10s} {real_d[v]:+11.3f} {a.mean():+10.3f} {m.mean():+13.3f} "
              f"{rep:10.0f}%")

    print(f"\n  additive check (dSSW - bias should recover all-events):")
    for v in sno:
        imp = real_d[v] - res["variables"][v]["null_rate_matched"]
        print(f"    {v:10s} implied {imp:+.3f}  vs all-events {real_all[v]:+.3f}  "
              f"residual {imp - real_all[v]:+.3f}")

    (RESULTS / "selection_snotel.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> selection_snotel.json")


if __name__ == "__main__":
    main()
