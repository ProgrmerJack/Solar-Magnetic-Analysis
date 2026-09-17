#!/usr/bin/env python3
"""
validate_R_diagnostic.py
========================
Turns R from an assertion into a calibrated instrument.

THE GAP THIS CLOSES
  causal_timescale_ratio.py reports R = within-winter / between-winter effect,
  and asserts R~1 means an event-scale (causal) association and R~0 means a
  winter-scale (common-cause) one. That follows from the CONSTRUCTION:

      between-winter effect = causal + winter_offset
      within-winter effect  = causal            (the winter offset cancels)
      => R = causal / (causal + winter_offset) = the causal fraction

  But that is the population identity. It says nothing about whether the
  ESTIMATOR recovers it from finite, autocorrelated, seasonally-varying data
  with day-of-year adjustment, 43 events, and a winter-block bootstrap. Reporting
  R = 0.82 as "82% causal" without checking that is exactly the kind of
  unvalidated inference this project has had to withdraw four times.

THE TEST
  Synthetic winters built to match the real AO:
    - AR(1) with phi = 0.943, the measured lag-1 autocorrelation of winter AO
    - an annual harmonic seasonal cycle
    - 77 winters, 43 events in 36 winters -- the observed configuration
  Then the association is split by construction into two parts whose relative
  size is SET, not estimated:

    CAUSAL       days +0..+60 after onset are depressed by c
    COMMON CAUSE winters that host an event are depressed by k THROUGHOUT,
                 which also makes those winters "SSW-prone" in the same way a
                 real confounder would

    true_f = c / (c + k)

  R is then estimated exactly as in causal_timescale_ratio.py and compared with
  true_f across the full range. What is measured:
    bias      does R track true_f, or is it systematically off?
    coverage  does the winter-block bootstrap interval contain true_f at 95%?
    power     can R distinguish f=0 from f=1 at n=43 events?

  If R is biased or its interval under-covers, the observational R = 0.82 cannot
  be read as a causal fraction and the claim must be weakened accordingly.

Output: validate_R_diagnostic.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
import causal_timescale_ratio as CTR                # noqa: E402

PHI = 0.943               # measured lag-1 autocorrelation of observed winter AO
N_WINTERS = 77
N_EVENT_WINTERS = 36
N_EVENTS = 43
WIN_LEN = 151             # 1 Nov .. 31 Mar
TOTAL = 1.0               # total association size, split between the two parts
N_REP = 300
F_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]


def make_winter(rng, n=WIN_LEN):
    """AR(1) series with the observed persistence, plus a seasonal cycle."""
    e = rng.normal(0, np.sqrt(1 - PHI ** 2), n)
    x = np.empty(n)
    x[0] = rng.normal(0, 1)
    for i in range(1, n):
        x[i] = PHI * x[i - 1] + e[i]
    doy = np.arange(n)
    x = x + 0.30 * np.sin(2 * np.pi * doy / 365.25) + 0.15 * np.cos(2 * np.pi * doy / 365.25)
    return x


def simulate(rng, true_f, where="post"):
    """Build one synthetic record with a known causal fraction.

    `where` places the event-scale anomaly:
      "post"  days 0..+60 after onset   -- a genuine downward response
      "pre"   days -60..0 before onset  -- a PRECURSOR, which cannot have been
              caused by the event that follows it
    The original version only ever simulated "post". That is exactly why R's
    causal reading survived: the design could not produce the one structure
    that would have falsified it.
    """
    c = TOTAL * true_f            # event-scale (causal) depression
    k = TOTAL * (1 - true_f)      # winter-scale (common-cause) depression
    start = pd.Timestamp("1950-11-01")
    rows = []
    ev_winters = set(rng.choice(np.arange(N_WINTERS), N_EVENT_WINTERS, replace=False))
    # distribute events across event-winters (a few winters get two)
    counts = {w: 1 for w in ev_winters}
    for w in rng.choice(sorted(ev_winters), N_EVENTS - N_EVENT_WINTERS, replace=False):
        counts[w] += 1
    onsets = []
    for w in range(N_WINTERS):
        t0 = start + pd.DateOffset(years=w)
        idx = pd.date_range(t0, periods=WIN_LEN, freq="D")
        y = make_winter(rng)
        if w in ev_winters:
            y = y - k                                  # whole winter depressed
            picks = sorted(rng.choice(np.arange(20, WIN_LEN - 70),
                                      counts[w], replace=False))
            # enforce >=20 d separation so bins stay well defined
            picks = [p for i, p in enumerate(picks)
                     if i == 0 or p - picks[i - 1] >= 25]
            for p in picks:
                if where == "post":
                    y[p:min(p + 61, WIN_LEN)] -= c     # post-onset depression
                else:
                    y[max(0, p - 60):p] -= c           # PRE-onset depression
                onsets.append(idx[p])
        rows.append(pd.DataFrame({"y": y, "doy": idx.dayofyear,
                                  "winter": 1950 + w + 1}, index=idx))
    d = pd.concat(rows)
    return d, pd.DatetimeIndex(sorted(onsets))


def run_grid(rng, res, WHERE):
    # R must be measured WHERE the anomaly is. A first version left the window
    # at the post-onset (15,29) for both arms, so the pre-onset arm measured a
    # window containing no anomaly: the between-winter denominator went to zero
    # and R became a garbage ratio (sd = 26 at f=1). That is a near-zero-
    # denominator artefact, not evidence that R separates pre from post.
    CTR.WINDOW = (15, 29) if WHERE == "post" else (-30, -16)
    print(f"  [measuring R at lag window {CTR.WINDOW}]")
    print(f"AR(1) phi={PHI}, {N_WINTERS} winters, {N_EVENTS} events in "
          f"{N_EVENT_WINTERS} winters -- the observed configuration")
    print(f"{N_REP} replicates per point\n")
    print(f"{'true f':>7s} {'mean R':>8s} {'bias':>8s} {'sd':>7s} "
          f"{'5-95 pct':>18s} {'coverage':>9s}")
    print("-" * 62)

    for f in F_GRID:
        ests, covs = [], []
        for _ in range(N_REP):
            d, on = simulate(rng, f, WHERE)
            pairs = CTR.per_event_pairs(d, on, ycol="y")
            if len(pairs) < 10:
                continue
            b = np.array([p[0] for p in pairs])
            w = np.array([p[1] for p in pairs])
            if abs(b.mean()) < 1e-6:
                continue
            R = w.mean() / b.mean()
            ests.append(R)
            # winter-block bootstrap interval for this replicate
            wid = np.array([p[2] for p in pairs])
            uw = np.unique(wid)
            rs = []
            for _ in range(200):
                pick = rng.choice(uw, len(uw), replace=True)
                sel = np.concatenate([np.flatnonzero(wid == x) for x in pick])
                bb, ww = b[sel].mean(), w[sel].mean()
                if abs(bb) > 1e-6:
                    rs.append(ww / bb)
            if rs:
                lo, hi = np.percentile(rs, [2.5, 97.5])
                covs.append(lo <= f <= hi)
        a = np.array(ests)
        cov = float(np.mean(covs)) if covs else np.nan
        res[str(f)] = {
            "mean_R": round(float(a.mean()), 4),
            "bias": round(float(a.mean() - f), 4),
            "sd": round(float(a.std()), 4),
            "pct5_95": [round(float(np.percentile(a, 5)), 3),
                        round(float(np.percentile(a, 95)), 3)],
            "coverage": round(cov, 3), "n_ok": int(len(a))}
        print(f"{f:7.2f} {a.mean():8.3f} {a.mean()-f:+8.3f} {a.std():7.3f} "
              f"[{np.percentile(a,5):+7.2f},{np.percentile(a,95):+7.2f}] {cov:9.3f}")

    b = [res[str(f)]["bias"] for f in F_GRID]
    c = [res[str(f)]["coverage"] for f in F_GRID]
    print(f"\n  max |bias| across the grid: {max(abs(x) for x in b):.3f}")
    print(f"  coverage range: {min(c):.3f}..{max(c):.3f}  (nominal 0.95)")
    ok = max(abs(x) for x in b) < 0.10 and min(c) > 0.90
    print(f"  -> R is {'a CALIBRATED estimator of the causal fraction' if ok else 'NOT adequately calibrated -- weaken the claim'}")

    # power: can f=0 be told from f=1 at n=43?
    a0, a1 = res["0.0"], res["1.0"]
    print()
    print(f"  separation at n={N_EVENTS}: f=0 gives R={a0['mean_R']:+.2f} "
          f"(5-95 {a0['pct5_95']}), f=1 gives R={a1['mean_R']:+.2f} "
          f"(5-95 {a1['pct5_95']})")
    print(f"  overlap: {'NO -- the diagnostic separates them' if a0['pct5_95'][1] < a1['pct5_95'][0] else 'YES -- limited power at this n'}")
    return res


def main():
    rng = np.random.default_rng(20260731)
    out = {"phi": PHI, "n_winters": N_WINTERS, "n_events": N_EVENTS,
           "n_rep": N_REP, "grid": F_GRID, "arms": {}}

    print("=" * 70)
    print("ARM 1 -- anomaly AFTER onset (the only structure ever simulated)")
    print("=" * 70)
    out["arms"]["post_onset"] = run_grid(rng, {}, "post")

    print()
    print("=" * 70)
    print("ARM 2 -- anomaly BEFORE onset: a PRECURSOR, true causal fraction 0")
    print("=" * 70)
    print("A precursor cannot be caused by the event that follows it. If R were")
    print("a causal fraction it would return ~0 here for every f. If R measures")
    print("temporal CONCENTRATION it will return the same ~0.8 as arm 1.")
    print()
    out["arms"]["pre_onset"] = run_grid(rng, {}, "pre")

    post = out["arms"]["post_onset"]["1.0"]["mean_R"]
    pre = out["arms"]["pre_onset"]["1.0"]["mean_R"]
    out["verdict"] = {
        "R_post_onset_at_f1": post,
        "R_pre_onset_at_f1": pre,
        "difference": round(pre - post, 4),
        "reading": (
            "R returns the same value for an anomaly placed BEFORE onset as for "
            "one placed AFTER it. A pre-onset anomaly has causal fraction zero "
            "by construction, so R is NOT an estimator of the causal fraction; "
            "it measures how tightly the anomaly is concentrated around the "
            "onset date. The original validation simulated only the post-onset "
            "structure and therefore could not detect this."
            if abs(pre - post) < 0.25 else
            "R distinguishes pre- from post-onset placement; the causal reading "
            "survives this test and the concentration interpretation is wrong.")}
    print()
    print("=" * 70)
    print(f"  R at f=1, anomaly AFTER  onset : {post:+.3f}")
    print(f"  R at f=1, anomaly BEFORE onset : {pre:+.3f}")
    print(f"  difference                     : {pre-post:+.3f}")
    print("=" * 70)
    print("  " + out["verdict"]["reading"])

    (RESULTS / "validate_R_diagnostic.json").write_text(
        json.dumps(out, indent=2), encoding="utf8")
    print()
    print("Saved -> validate_R_diagnostic.json")


if __name__ == "__main__":
    main()
