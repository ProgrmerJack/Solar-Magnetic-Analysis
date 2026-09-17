#!/usr/bin/env python3
"""
granger_direction.py
====================
Tests the DIRECTION of stratosphere-surface coupling, which R cannot.

WHY THIS IS NEEDED
  R = within-winter/between-winter is flat across lag: 0.79 before onset, 0.82
  after. A precursor cannot be caused by the event that follows it, so R~0.8 at
  NEGATIVE lag proves R is not a causal fraction -- it measures how temporally
  CONCENTRATED the anomaly is around the event. That still rejects a winter-scale
  common cause (p < 0.0007), which is real, but a short-timescale process driving
  both the vortex and the surface would produce the same signature.

  Causation, unlike correlation or concentration, is ASYMMETRIC IN TIME. That is
  the property R throws away and this exploits.

THE TEST
  Granger causality on daily winter series, both directions:

      does u10(t-1..t-L) improve prediction of AM(t) given AM's own past?   (down)
      does AM(t-1..t-L) improve prediction of u10(t) given u10's own past?   (up)

  Downward coupling predicts a strong down->up asymmetry. A shared fast driver
  with no directional influence predicts symmetry.

  CONDITIONAL version: the same, additionally controlling for the tropospheric
  wave driving v'T' at 100 hPa and its lags. Wave activity drives the vortex AND
  projects on the surface, so it is the obvious common cause; conditioning on it
  is what separates "stratosphere influences surface" from "wave activity does
  both". This is the test the mediation analysis approaches from the other side.

  Inference is by winter-block bootstrap, never by the F-distribution: daily
  series are strongly autocorrelated (AO lag-1 = 0.943) and asymptotic p-values
  would be badly anti-conservative.

Output: granger_direction.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "5_mechanism"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import canonical_event_study as K                   # noqa: E402
import ensemble_precursor as EP                     # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
LAGS = 10
# Raised from 800 to 2000 on 2026-09-17. The project's standing floor is 1000
# (FAILURES.md, 2026-07-28) and this result was the only live one below it;
# MPI-ESM-1-2-HAM sat at P=0.063, close enough to 0.05 to need the headroom.
N_BOOT = 2000
SEASON = (11, 12, 1, 2, 3)


def lagmat(x, L):
    return np.column_stack([np.roll(x, i) for i in range(1, L + 1)])


def gc(y, x, winter, z=None, L=LAGS):
    """Variance ratio for x -> y, given y's own past (and optionally z's)."""
    Y = y[L:]
    base = [lagmat(y, L)[L:]]
    if z is not None:
        base.append(lagmat(z, L)[L:])
    full = base + [lagmat(x, L)[L:]]
    w = winter[L:]
    # drop rows spanning a winter boundary
    ok = np.ones(len(Y), bool)
    for i in range(1, L + 1):
        ok &= (np.roll(winter, i)[L:] == w)
    Y, w = Y[ok], w[ok]
    A0 = np.column_stack([np.ones(ok.sum())] + [b[ok] for b in base])
    A1 = np.column_stack([np.ones(ok.sum())] + [b[ok] for b in full])

    def rss(A, yy):
        b, *_ = np.linalg.lstsq(A, yy, rcond=None)
        r = yy - A @ b
        return float(r @ r)
    return rss(A0, Y), rss(A1, Y), w, A0, A1, Y


def stat(y, x, winter, z=None):
    r0, r1, w, A0, A1, Y = gc(y, x, winter, z)
    return 1.0 - r1 / r0, (w, A0, A1, Y)


def boot_diff(y, x, winter, z=None, seed=1):
    """Bootstrap the DOWN-minus-UP asymmetry over winters."""
    sd, (w, A0d, A1d, Yd) = stat(y, x, winter, z)
    su, (w2, A0u, A1u, Yu) = stat(x, y, winter, z)
    rng = np.random.default_rng(seed)
    uw = np.unique(w)
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(uw, len(uw), replace=True)
        sel = np.concatenate([np.flatnonzero(w == q) for q in pick])
        sel2 = np.concatenate([np.flatnonzero(w2 == q) for q in pick])

        def r2(A0, A1, Y, s):
            b0, *_ = np.linalg.lstsq(A0[s], Y[s], rcond=None)
            b1, *_ = np.linalg.lstsq(A1[s], Y[s], rcond=None)
            e0 = Y[s] - A0[s] @ b0
            e1 = Y[s] - A1[s] @ b1
            return 1.0 - float(e1 @ e1) / float(e0 @ e0)
        out.append(r2(A0d, A1d, Yd, sel) - r2(A0u, A1u, Yu, sel2))
    o = np.array(out)
    return sd, su, o


def obs_series():
    ao = K.load_series("AO")["y"]
    st = pd.read_parquet(STRAT)
    if st.index.tz is not None:
        st.index = st.index.tz_convert("UTC").tz_localize(None)
    u = st["uwnd_ms_10hPa"]
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"]
    idx = ao.index.intersection(u.index).intersection(hf.index)
    idx = idx[np.isin(idx.month, SEASON)]
    d = pd.DataFrame({"am": ao.reindex(idx), "u10": u.reindex(idx),
                      "vt": hf.reindex(idx)}, index=idx).dropna()
    d["winter"] = EP.winter_of(d.index)
    for c in ("am", "u10", "vt"):
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    return d


def report(name, d, use_z):
    z = d["vt"].values if use_z else None
    sd, su, o = boot_diff(d["am"].values, d["u10"].values, d["winter"].values, z)
    p = float((o <= 0).mean())
    lab = "conditional on wave driving" if use_z else "bivariate"
    print(f"  {name:16s} {lab:28s} down {sd:.4f}  up {su:.4f}  "
          f"asym {sd-su:+.4f} [{np.percentile(o,2.5):+.4f},"
          f"{np.percentile(o,97.5):+.4f}]  P(asym<=0)={p:.4f}")
    return {"down": round(sd, 5), "up": round(su, 5),
            "asymmetry": round(sd - su, 5),
            "asym_CI95": [round(float(np.percentile(o, 2.5)), 5),
                          round(float(np.percentile(o, 97.5)), 5)],
            "p_asym_le_0": p, "conditional": use_z}


def main():
    out = {"lags": LAGS, "n_boot": N_BOOT, "results": {}}
    print("Granger asymmetry: (vortex -> surface) minus (surface -> vortex)")
    print("  positive asymmetry = downward influence dominates\n")
    d = obs_series()
    print(f"  observations: {len(d):,} winter days, {d['winter'].nunique()} winters")
    out["results"]["observations_bivariate"] = report("OBSERVATIONS", d, False)
    out["results"]["observations_conditional"] = report("OBSERVATIONS", d, True)

    for name in ("CanESM5", "IPSL-CM5A2-INCA", "MPI-ESM-1-2-HAM", "CESM2-FV2"):
        files = sorted(RAW.glob(f"{name}_*_zm.nc"))[:3]
        if not files:
            continue
        fr = []
        for k, f in enumerate(files):
            m = EP.load_member(f)
            m = m[np.isin(m.index.month, SEASON)].dropna()
            m = m.assign(winter=EP.winter_of(m.index) + 100000 * (k + 1))
            fr.append(m)
        dm = pd.concat(fr)
        for c in ("am", "u10"):
            dm[c] = (dm[c] - dm[c].mean()) / dm[c].std()
        out["results"][f"{name}_bivariate"] = report(name, dm, False)

    o = out["results"]["observations_conditional"]
    print(f"\n=== VERDICT ===")
    print(f"  conditional on wave driving, the vortex->surface direction is "
          f"{'STRONGER' if o['asymmetry']>0 else 'NOT stronger'} than the reverse")
    print(f"  asymmetry {o['asymmetry']:+.4f}, P(<=0) = {o['p_asym_le_0']:.4f}")
    print("  a shared fast driver with no directional influence predicts asym = 0")
    (RESULTS / "granger_direction.json").write_text(json.dumps(out, indent=2),
                                                 encoding="utf8")
    print("\nSaved -> granger_direction.json")


if __name__ == "__main__":
    main()
