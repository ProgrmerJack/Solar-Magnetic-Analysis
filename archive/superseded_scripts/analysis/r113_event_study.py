#!/usr/bin/env python3
"""
r113_event_study.py
===================
Rebuild of R112 as a PROPER EVENT STUDY, after peer review identified a genuine
identification flaw.

THE FLAW IN R112
  R112 fitted one lag-window indicator at a time:
        y_d = alpha_w + beta*W_d + harmonics + eps
  so the comparison group for ANY window was "all other days in that winter".
  When estimating the -45..-16 placebo, the controls therefore included the
  central date and the whole 0..+60 response period. A null placebo coefficient
  did not establish "the estimator reports nothing before onset"; it only showed
  that the pre-onset window did not differ from a control group that was itself
  contaminated by the response. The same criticism applies in reverse to the
  post-onset estimate.

THE FIX
  Estimate mutually exclusive lead and lag bins SIMULTANEOUSLY against an
  explicit omitted baseline of winter days more than `BASELINE_GAP` days from
  ANY SSW central date:

        y_d = alpha_w + sum_b beta_b * 1[d in bin b] + harmonics + eps

  Days are assigned to the bin of the NEAREST central date, so bins are mutually
  exclusive even when a winter contains two SSWs (as 2023/24 does). Days falling
  in no bin and outside the washout are the baseline.

  With this specification beta_b for a pre-onset bin is a genuine placebo: it
  compares pre-onset days with days far from any event, not with the response.

ADDITIONAL VALIDATION REQUESTED BY REVIEW
  * pseudo-onset randomisation: re-draw fake central dates within the same
    winters and refit, giving a null distribution for the whole profile;
  * effect injection: add a known synthetic response to the observed series and
    confirm the estimator recovers its magnitude (guards against the positive
    control being circular, since the Baldwin-Dunkerton magnitude itself comes
    from composites);
  * baseline-gap sensitivity.

Output: data/results/r113_event_study.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r113_event_study.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_HARM = 3
N_BOOT = 2000
N_PSEUDO = 1000
BASELINE_GAP = 75        # days from any onset to count as baseline

BINS = [(-60, -46), (-45, -31), (-30, -16), (-15, -1),
        (0, 14), (15, 29), (30, 44), (45, 60)]
BIN_LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]


def read_cpc(path):
    d = pd.read_csv(path, sep=r"\s+", header=None,
                    names=["year", "month", "day", "value"])
    d["date"] = pd.to_datetime(dict(year=d.year, month=d.month, day=d.day),
                               errors="coerce")
    return d.dropna(subset=["date"]).set_index("date")["value"].astype(float).sort_index()


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def build(s, onsets):
    """Daily frame with nearest-event lag and mutually exclusive bin codes."""
    d = pd.DataFrame({"y": s.values}, index=s.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    doy = d.index.dayofyear.values
    for k in range(1, N_HARM + 1):
        d[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        d[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)

    on = np.sort(np.array([np.datetime64(o, "D") for o in onsets]))
    dd = d.index.values.astype("datetime64[D]")
    # signed lag to the NEAREST onset -> guarantees mutually exclusive bins
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    d["lag"] = lag

    code = np.full(len(d), -1)          # -1 = not in any bin
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    d["bin"] = code
    d["is_baseline"] = np.abs(lag) > BASELINE_GAP
    # drop washout days: outside all bins but within BASELINE_GAP of an event
    d = d[(d["bin"] >= 0) | d["is_baseline"]].copy()
    return d


def design(d):
    X = np.column_stack(
        [(d["bin"].values == i).astype(float) for i in range(len(BINS))]
        + [d[f"{t}{k}"].values for k in range(1, N_HARM + 1) for t in ("s", "c")])
    return d["y"].values.astype(float), X, pd.Categorical(d["winter"]).codes.astype(np.int64)


def _demean(A, codes, n):
    cnt = np.bincount(codes, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        sm = np.bincount(codes, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (sm / cnt)[codes]
    return out


def fit(y, X, codes):
    _, c = np.unique(codes, return_inverse=True)
    A = np.column_stack([y, X])
    Ad = _demean(A, c, c.max() + 1)
    try:
        beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return beta[:len(BINS)]
    except Exception:
        return None


def run(d, n_boot=N_BOOT, seed=0):
    y, X, codes = design(d)
    est = fit(y, X, codes)
    if est is None:
        return None
    rng = np.random.default_rng(seed)
    _, c = np.unique(codes, return_inverse=True)
    n = c.max() + 1
    idx_by_w = [np.flatnonzero(c == w) for w in range(n)]
    boot = []
    for _ in range(n_boot):
        pick = rng.integers(0, n, n)
        rows = np.concatenate([idx_by_w[k] for k in pick])
        cb = np.concatenate([np.full(len(idx_by_w[k]), j) for j, k in enumerate(pick)])
        v = fit(y[rows], X[rows], cb)
        if v is not None and np.all(np.isfinite(v)):
            boot.append(v)
    boot = np.array(boot)
    out = {}
    for i, lab in enumerate(BIN_LABELS):
        b = boot[:, i]
        out[lab] = {
            "effect": round(float(est[i]), 4),
            "CI95": [round(float(np.percentile(b, 2.5)), 4),
                     round(float(np.percentile(b, 97.5)), 4)],
            "p_two_sided": float(2 * min((b >= 0).mean(), (b <= 0).mean())),
        }
    return out, est


def pseudo_onset_null(s, onsets, n=N_PSEUDO, seed=1):
    """Re-draw fake onsets inside the SAME winters; refit; null distribution."""
    rng = np.random.default_rng(seed)
    real_w = pd.DatetimeIndex(onsets)
    winters = winter_of(real_w)
    peak = np.zeros(n)
    for t in range(n):
        fake = []
        for w in winters:
            # a random mid-winter date in the same winter
            start = pd.Timestamp(year=int(w) - 1, month=12, day=1)
            fake.append(start + pd.Timedelta(days=int(rng.integers(0, 100))))
        d = build(s, pd.DatetimeIndex(fake))
        y, X, codes = design(d)
        b = fit(y, X, codes)
        if b is not None:
            peak[t] = b[4:7].mean()      # mean of the 0..+44 post-onset bins
    return peak


def injection_test(s, onsets, amplitude=-0.5, seed=2):
    """Add a known synthetic post-onset response; check it is recovered."""
    s2 = s.copy()
    idx = s2.index.values.astype("datetime64[D]")
    add = np.zeros(len(s2))
    for o in onsets:
        o64 = np.datetime64(pd.Timestamp(o), "D")
        lag = (idx - o64).astype(int)
        add[(lag >= 0) & (lag <= 60)] += amplitude
    s2.iloc[:] = s2.values + add
    d = build(s2, onsets)
    y, X, codes = design(d)
    b = fit(y, X, codes)
    post = float(np.mean(b[4:7])) if b is not None else np.nan
    return {"injected_amplitude": amplitude,
            "recovered_post_mean": round(post, 4),
            "recovery_ratio": round(post / amplitude, 3) if amplitude else None}


def main():
    global BASELINE_GAP
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    onsets = pd.DatetimeIndex(pd.to_datetime(can["date"]).dropna())
    ao = read_cpc(ROOT / "data/processed/atmospheric/ao_daily_cpc.txt")
    nao = read_cpc(ROOT / "data/processed/atmospheric/nao_daily_cpc.txt")

    res = {"bins": BIN_LABELS, "baseline_gap_days": BASELINE_GAP,
           "n_onsets": len(onsets), "series": {}}

    for name, s in (("AO", ao), ("NAO", nao)):
        ev = onsets[(onsets >= s.index.min()) & (onsets <= s.index.max())]
        d = build(s, ev)
        nb = int(d["is_baseline"].sum())
        print(f"\n=== {name}: {len(ev)} events, {d['winter'].nunique()} winters, "
              f"{len(d):,} days ({nb:,} baseline) ===")
        r, est = run(d)
        for lab in BIN_LABELS:
            v = r[lab]
            flag = " *" if v["p_two_sided"] < 0.05 else ""
            print(f"  {lab:>10s}  {v['effect']:+7.3f}  "
                  f"[{v['CI95'][0]:+.3f},{v['CI95'][1]:+.3f}]  P={v['p_two_sided']:.4f}{flag}")
        res["series"][name] = {"n_events": len(ev),
                               "n_winters": int(d["winter"].nunique()),
                               "n_days": int(len(d)), "n_baseline_days": nb,
                               "bins": r}

    # validation on AO
    ev = onsets[(onsets >= ao.index.min()) & (onsets <= ao.index.max())]
    inj = injection_test(ao, ev)
    print(f"\ninjection test: injected {inj['injected_amplitude']}, "
          f"recovered {inj['recovered_post_mean']} "
          f"(ratio {inj['recovery_ratio']})")
    res["injection_test_AO"] = inj

    print("pseudo-onset randomisation (1000 draws) ...")
    null = pseudo_onset_null(ao, ev)
    d = build(ao, ev)
    y, X, codes = design(d)
    obs = float(np.mean(fit(y, X, codes)[4:7]))
    p = float((null <= obs).mean())
    print(f"  observed post-onset mean {obs:+.3f}; "
          f"pseudo-onset null mean {null.mean():+.3f} sd {null.std():.3f}; "
          f"P={p:.4f}")
    res["pseudo_onset_AO"] = {"observed_post_mean": round(obs, 4),
                              "null_mean": round(float(null.mean()), 4),
                              "null_sd": round(float(null.std()), 4),
                              "p_one_sided": p, "n_draws": int(len(null))}

    # baseline-gap sensitivity
    sens = {}
    for gap in (60, 75, 90):
        BASELINE_GAP = gap
        d = build(ao, ev)
        y, X, codes = design(d)
        b = fit(y, X, codes)
        sens[gap] = {"pre_-45..-31": round(float(b[1]), 4),
                     "post_0..+14": round(float(b[4]), 4),
                     "post_+15..+29": round(float(b[5]), 4)}
    BASELINE_GAP = 75
    res["baseline_gap_sensitivity_AO"] = sens
    print("baseline-gap sensitivity:", json.dumps(sens))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
