#!/usr/bin/env python3
"""
enso_confound.py
================
The AAO negative control fails and the estimator is not to blame:

  observed      -0.659 at +45..+60 d, bootstrap p=0.0147
  seed-stable   significant on 60/60 fixed seeds at 6,000 replicates
  estimator OK  clean pseudo-onset null centred on zero (mean -0.018)
  bootstrap OK  like-for-like calibration at that bin: SE ratio 1.06,
                coverage 0.945 (boot_calibration.py)

An earlier reading of this file claimed the control PASSED because a
randomisation p of 0.081 disagreed with the bootstrap p of 0.0147. That was
wrong: the null SD came from 5,121 clean days and the bootstrap from 8,638, and
sqrt(8638/5121)=1.30 rescales 0.368 to 0.283 against a bootstrap SE of 0.268.
The two agree. The gap was sample size, exactly as in the VOID ao_null_perbin.py.

So something genuinely links NH SSW onsets to the Southern annular mode.

HYPOTHESIS: ENSO is a common cause. It modulates NH polar vortex variability
(and so SSW likelihood) and independently projects onto the SAM. If that is the
mechanism, the AAO is not a valid negative control for SSW event studies and
Gate 6's PREMISE is wrong, rather than the estimator.

TEST: add ONI as a covariate to the event-study regression, interpolated to
daily and entered alongside the harmonics and bin dummies. Compare the +45..+60
AAO effect with and without it. The AO is run through the same treatment as a
control -- if adding ONI also guts the AO response, the covariate is absorbing
signal rather than confounding, and the test is uninformative.

Output: enso_confound.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "2_event_study"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue
import calibrate_seasonality as C
import multi_index_event_study as M

ONI = HERE.parents[0] / "03_data_ingestion" / "enso_oni.parquet"
N_BOOT = 6000
WATCH = {"aao": "+45..+60", "ao": "+15..+29"}


def load_oni(index):
    o = pd.read_parquet(ONI).set_index("date")["oni"].sort_index()
    # monthly -> daily by time interpolation; ONI is already a 3-month running mean
    daily = o.reindex(o.index.union(index)).interpolate("time").reindex(index)
    return daily.values.astype(float)


def fit_with(dd, oni):
    y = dd["y"].values.astype(float)
    S = C.harmonics(dd["doy"].values, 3)
    D = np.column_stack([(dd["bin"].values == i).astype(float) for i in range(len(M.BINS))])
    cols = [y, D, S] if oni is None else [y, D, S, oni[:, None]]
    A = M.demean(np.column_stack(cols), dd["winter"].values)
    b, *_ = np.linalg.lstsq(A[:, 1:], A[:, 0], rcond=None)
    return b[:len(M.BINS)], (None if oni is None else float(b[-1]))


def run(name):
    d = M.load(name)
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    dd = M.build(d, on)
    oni_all = load_oni(dd.index)
    ok = np.isfinite(oni_all)
    if not ok.all():
        dd, oni_all = dd[ok], oni_all[ok]
    out = {"n_events": int(len(on)), "watch_bin": WATCH[name], "variants": {}}
    print(f"\n=== {name.upper()}  n={len(on)} events, {len(dd)} days, "
          f"ONI sd {np.std(oni_all):.2f} ===")

    for lab, oni in (("without_ONI", None), ("with_ONI", oni_all)):
        est, coef = fit_with(dd, oni)
        rng = np.random.default_rng(20260730)
        winters = np.array(sorted(dd["winter"].unique()))
        idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
        boot = []
        for _ in range(N_BOOT):
            pick = rng.choice(winters, len(winters), replace=True)
            rows = np.concatenate([idx[w] for w in pick])
            sub = dd.iloc[rows].copy()
            sub["winter"] = np.concatenate(
                [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
            try:
                b, _ = fit_with(sub, None if oni is None else oni[rows])
                if np.all(np.isfinite(b)):
                    boot.append(b)
            except Exception:
                pass
        boot = np.array(boot)
        j = M.LABELS.index(WATCH[name])
        col = boot[:, j]
        p = float(2 * min((col >= 0).mean(), (col <= 0).mean()))
        out["variants"][lab] = {
            "effect": round(float(est[j]), 4),
            "CI95": [round(float(np.percentile(col, 2.5)), 4),
                     round(float(np.percentile(col, 97.5)), 4)],
            "p_two_sided": p, "oni_coefficient": None if coef is None else round(coef, 4)}
        print(f"  {lab:12s} {WATCH[name]} = {est[j]:+.3f} "
              f"[{np.percentile(col,2.5):+.3f},{np.percentile(col,97.5):+.3f}] "
              f"p={p:.4f}" + ("" if coef is None else f"   ONI beta {coef:+.4f}"))
    a = out["variants"]["without_ONI"]["effect"]
    b = out["variants"]["with_ONI"]["effect"]
    out["attenuation_fraction"] = round(float(1 - b / a), 3) if a else None
    print(f"  -> effect changes {a:+.3f} -> {b:+.3f} "
          f"({100*(1-b/a):.0f}% attenuation)")
    return out


def main():
    res = {"n_boot": N_BOOT, "series": {n: run(n) for n in ("aao", "ao")}}
    aao = res["series"]["aao"]; ao = res["series"]["ao"]
    print("\n=== verdict ===")
    print(f"  AAO attenuation {100*aao['attenuation_fraction']:.0f}%, "
          f"AO attenuation {100*ao['attenuation_fraction']:.0f}%")
    print("  AAO killed while AO survives -> ENSO is a common cause; the AAO is")
    print("     NOT a valid negative control and Gate 6's premise is wrong")
    print("  both killed -> the covariate absorbs signal; test uninformative")
    print("  neither killed -> ENSO is not the mechanism; failure unexplained")
    (RESULTS / "enso_confound.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> enso_confound.json")


if __name__ == "__main__":
    main()
