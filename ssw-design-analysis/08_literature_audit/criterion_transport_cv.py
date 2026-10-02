#!/usr/bin/env python3
"""
criterion_transport_cv.py
=========================
DOES A SHIFT ESTIMATED WITHOUT AN EVENT REPRODUCE THAT EVENT'S OUTCOME AND LABEL?
WHOLE-WINTER CROSS-VALIDATION OF THE SHIFTED-EVENT-FREE NULL (ERA5).

Plan approved 2026-10-02 (user chose to act on the external audit, section 5).
Committed BEFORE it was run. RETROSPECTIVE cross-validation: the 39 outcomes have
been examined before (criterion_sweep.py, era5_regional_test.py); this is not new
untouched evidence, only a test that the null transports when its shift is
estimated without the evaluated winter.

DESIGN
  Events: the 39 catalogued SSWs with ERA5 NAM outcomes (criterion_sweep's set).
  Outer folds: hindcast winters; every event of the held-out winter is predicted
  from the other winters only. In each fold, from the training winters alone:
    event-free candidates (gate3 zone-free days, training winters only);
    the shift of each level = mean window-mean NAM of the training SSWs minus the
      mean over 200 day-of-year-matched event-free sets drawn from training
      winters (criterion_sweep's estimator, restricted to training data).
  Forecasts for a held-out event (day of year d): the training-winter event-free
  days within +-10 days of d, each displaced by the training shift, form an
  ensemble; its label probability is the share classified downward and its
  window means the predictive distribution.
  Competitors on the same events: unshifted climatology (the same ensemble, no
  shift) and a constant class rate (the training SSWs' downward share).
  Primary criterion: the published surface conditions (1000 hPa NAM, days 8-52,
  window mean < 0 and > 50% of days < 0); secondary: with condition 3 (150 hPa);
  and all 108 versions of criterion_sweep summarised.
STATISTICS
  Calibration: observed downward count against the sum of predicted
  probabilities; Poisson-binomial two-sided p (100,000 simulations).
  Brier score of the shifted null against climatology and against the constant
  rate; CRPS (ensemble form) of the window mean, shifted against unshifted; mean
  PIT and share in the outer deciles. Score differences with 2,000 bootstrap
  resamples of winters.
READING (fixed now): the null transports if the calibration p > 0.05 and the
  shifted null's CRPS beats climatology (interval excluding zero); it adds nothing
  over a constant class rate if their Brier scores do not differ (interval
  including zero) -- which is what one shifted population predicts for the label.

Output: results/current/9_literature/criterion_transport_cv.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import criterion_sweep as CS                         # noqa: E402

G, load_catalogue = CS.G, CS.load_catalogue
NAME = "criterion_transport_cv"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RESULTS = CS.RESULTS
N_SHIFT_SETS, N_BOOT, N_SIM, HALFWIN = 200, 2000, 100000, 10
PRIMARY = (52, "nam_1000", 0.0, 0.5, False)
SECONDARY = (52, "nam_1000", 0.0, 0.5, True)


def winter(t):
    return t.year if t.month <= 6 else t.year + 1


def crps_ens(x, y):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.mean(np.abs(x - y)) - 0.5 * np.mean(np.abs(x[:, None] - x[None, :])))


def main():
    rng = np.random.default_rng(SEED)
    e = pd.read_parquet(CS.ERA5)
    allev = load_catalogue("primary")
    real = allev[(allev >= e.index.min() + pd.Timedelta(days=70)) & (allev <= e.index.max() - pd.Timedelta(days=70))]
    clean = G.zone_free_index(e.index, allev)
    levels = (*CS.LEVELS, "nam_150")
    Mr = {c: CS.window_matrix(e[c], real) for c in levels}
    Mc = {c: CS.window_matrix(e[c], clean) for c in levels}
    wr = np.array([winter(t) for t in real]); wc = np.array([winter(t) for t in clean])
    dr = np.array([t.dayofyear for t in real]); dc = np.array([t.dayofyear for t in clean])
    specs = CS.specs(CS.LEVELS, (False, True))
    folds = sorted(set(wr))
    pred = {s: np.full(len(real), np.nan) for s in specs}
    rate = {s: np.full(len(real), np.nan) for s in specs}
    unsh = {s: np.full(len(real), np.nan) for s in specs}
    obs = {s: np.full(len(real), np.nan) for s in specs}
    crps_s = np.full(len(real), np.nan); crps_c = np.full(len(real), np.nan); pit = np.full(len(real), np.nan)
    for w in folds:
        te = np.flatnonzero(wr == w); tr = np.flatnonzero(wr != w)
        ctr = np.flatnonzero(wc != w)
        clean_tr = clean[ctr]; by = G.build_by_doy(clean_tr)
        pos = {t: i for i, t in enumerate(clean)}
        shift = {}
        sets = [G.draw_clean(clean_tr, dr[tr], rng, by) for _ in range(N_SHIFT_SETS)]
        for W in CS.WS:
            for c in levels:
                pool = np.concatenate([CS.wmean(Mc[c][[pos[t] for t in s_]], W) for s_ in sets])
                shift[(W, c)] = float(np.nanmean(CS.wmean(Mr[c][tr], W)) - np.nanmean(pool))
        for i in te:
            dd = np.abs(dc[ctr] - dr[i]); dd = np.minimum(dd, 366 - dd)
            cand = ctr[dd <= HALFWIN]
            for s in specs:
                W, L, tau, f, c3 = s
                lab_r, ok_r = CS.classify(Mr[L][[i]], Mr["nam_150"][[i]], W, tau, f, c3)
                if not ok_r[0]:
                    continue
                lab_c, ok_c = CS.classify(Mc[L][cand], Mc["nam_150"][cand], W, tau, f, c3,
                                          shift[(W, L)], shift[(W, "nam_150")] if c3 else 0.0)
                lt, okt = CS.classify(Mr[L][tr], Mr["nam_150"][tr], W, tau, f, c3)
                lab_u, ok_u = CS.classify(Mc[L][cand], Mc["nam_150"][cand], W, tau, f, c3)
                pred[s][i] = lab_c[ok_c].mean(); rate[s][i] = lt[okt].mean(); obs[s][i] = float(lab_r[0])
                unsh[s][i] = lab_u[ok_u].mean()
            y = CS.wmean(Mr["nam_1000"][[i]], 52)[0]
            ens = CS.wmean(Mc["nam_1000"][cand], 52)
            if np.isfinite(y):
                crps_s[i] = crps_ens(ens + shift[(52, "nam_1000")], y); crps_c[i] = crps_ens(ens, y)
                xs = ens[np.isfinite(ens)] + shift[(52, "nam_1000")]
                pit[i] = (np.sum(xs < y) + 0.5) / (len(xs) + 1)
    res = {"plan_approved": "2026-10-02", "seed": SEED, "n_boot": N_BOOT, "n_events": int(len(real)),
           "n_winters": len(folds), "halfwin_days": HALFWIN, "label": "retrospective cross-validation (outcomes seen before)"}

    def evaluate(s):
        ok = np.isfinite(obs[s]) & np.isfinite(pred[s])
        p, r_, y, wv, uc = pred[s][ok], rate[s][ok], obs[s][ok], wr[ok], unsh[s][ok]
        sim = (rng.random((N_SIM, len(p))) < p).sum(1)
        k = int(y.sum())
        pcal = float(min(1, 2 * min(np.mean(sim <= k), np.mean(sim >= k))))
        uw = np.unique(wv)
        def sc(idx):
            return (np.mean((p[idx] - y[idx]) ** 2), np.mean((r_[idx] - y[idx]) ** 2))
        b = []
        for _ in range(N_BOOT):
            idx = np.concatenate([np.flatnonzero(wv == u) for u in rng.choice(uw, len(uw))])
            a, c = sc(idx); b.append(c - a)
        bs, br = sc(np.arange(len(y)))
        # registered but omitted in the first run (added 2026-10-02 after the code review):
        # Brier score of the shifted null against unshifted climatology
        bu = float(np.mean((uc - y) ** 2))
        bb = [np.mean((uc[ix] - y[ix]) ** 2) - np.mean((p[ix] - y[ix]) ** 2)
              for ix in (np.concatenate([np.flatnonzero(wv == u_) for u_ in rng.choice(uw, len(uw))]) for _ in range(N_BOOT))]
        return {"n": int(ok.sum()), "observed_dw": k, "expected_dw": round(float(p.sum()), 2),
                "p_calibration": round(pcal, 4), "brier_shifted": round(float(bs), 4), "brier_constant_rate": round(float(br), 4),
                "brier_rate_minus_shifted": round(float(br - bs), 4), "ci95": [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])],
                "mean_predicted": round(float(p.mean()), 4), "constant_rate_mean": round(float(r_.mean()), 4),
                "brier_climatology_unshifted": round(bu, 4), "brier_climatology_minus_shifted": round(float(bu - bs), 4),
                "ci95_climatology_minus_shifted": [round(float(x), 4) for x in np.percentile(bb, [2.5, 97.5])]}

    res["primary_published_surface"] = evaluate(PRIMARY)
    res["secondary_with_150hPa"] = evaluate(SECONDARY)
    allv = [evaluate(s) for s in specs]
    res["all_108_versions"] = {"n_versions": len(allv), "share_p_calibration_lt_0.05": round(float(np.mean([a["p_calibration"] < 0.05 for a in allv])), 4),
                               "median_expected_minus_observed": round(float(np.median([a["expected_dw"] - a["observed_dw"] for a in allv])), 2)}
    ok = np.isfinite(crps_s) & np.isfinite(crps_c)
    d = crps_c[ok] - crps_s[ok]; wv = wr[ok]; uw = np.unique(wv)
    bd = [np.mean(np.concatenate([d[wv == u] for u in rng.choice(uw, len(uw))])) for _ in range(N_BOOT)]
    res["crps_window_mean_days8_52"] = {"n": int(ok.sum()), "crps_shifted": round(float(crps_s[ok].mean()), 4),
                                        "crps_climatology": round(float(crps_c[ok].mean()), 4),
                                        "climatology_minus_shifted": round(float(d.mean()), 4),
                                        "ci95": [round(float(x), 4) for x in np.percentile(bd, [2.5, 97.5])],
                                        "mean_pit": round(float(np.nanmean(pit)), 4),
                                        "share_pit_outer_deciles": round(float(np.mean((pit[ok] < 0.1) | (pit[ok] > 0.9))), 4)}
    for k_ in ("primary_published_surface", "secondary_with_150hPa", "all_108_versions", "crps_window_mean_days8_52"):
        print(k_, json.dumps(res[k_]), flush=True)
    (RESULTS / "criterion_transport_cv.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> criterion_transport_cv.json")


if __name__ == "__main__":
    main()
