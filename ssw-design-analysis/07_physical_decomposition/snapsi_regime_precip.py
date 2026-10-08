#!/usr/bin/env python3
"""
snapsi_regime_precip.py
=======================
DOES A CHANGE IN HOW OFTEN CIRCULATION REGIMES OCCUR EXPLAIN THE SSW'S EFFECT ON
RAINFALL, OR MUST RAINFALL WITHIN THE REGIMES CHANGE?

Registered 2026-10-08, committed BEFORE any SNAPSI precipitation or regional pressure
value was fetched (03_data_ingestion/acquire_snapsi_regime_precip.py). This is the
physical alternative a class reading implies: the SSW acts by changing the frequency of
fixed weather types. It is tested here with types defined from circulation alone.

AMENDMENT 2026-10-08, before any rainfall outcome was computed (the fetch had read values
  only for range checks): the 2019 initialisations cover post-onset days 1-24 (s20181213)
  and 6-25 (s20190108) only -- the 8 January run starts after the 2 January onset -- so
  "days 1-25" below becomes DAYS 8-24 throughout (the paper's regional-temperature window,
  covered by all four initialisations), including the anomaly reference and every outcome.
  NRL is excluded before any outcome: its pr files state kg m-2 s-1 but their global mean
  converts to ~4,600 mm/day and varies erratically between steps, so the units cannot be
  determined (the fetch's plausibility check rejects every NRL member).

DATA  SNAPSI nudged and control, the four Northern Hemisphere initialisations
  (s20180125, s20180208 -> 2018 event; s20181213, s20190108 -> 2019 event), every centre
  with both arms; per member, days 1-25 after the event's onset (lead offset as
  snapsi_law.post_onset_offset):
    pressure  daily-mean psl (6hrPt), coarse-grained to 5 deg over 30-80N, 80W-40E;
    rainfall  daily-mean pr (6hr) over the boxes of Dai et al. 2025: IBERIA 36-44N,
              10-1W; SCANDINAVIA 55-70N, 4-20E; BRITISH_ISLES 49-61N, 10W-1E (mm/day).
REGIMES (fixed now; computed from pressure only, before any rainfall is read)
  Per centre: anomalies from that pair's control ensemble-mean field over days 1-25;
  EOFs (cos-latitude weighted) of the control anomalies pooled over its initialisations;
  the first 10 principal components; k-means, k = 4 (the Euro-Atlantic convention),
  50 restarts, seed fixed; nudged days assigned to the nearest control centroid in the
  same space. Sensitivity: k = 6.
FREQUENCY-ONLY CONSTRUCTION  For every nudged member and day d with regime j, rainfall
  (all boxes jointly) of a random control day of the same pair with regime j and day
  within d +- 5 (widened to +- 10, then the whole window, if none); the nudged member's
  regime sequence -- and so regime persistence -- is kept, the control's rainfall
  within regimes is kept. 200 synthetic ensembles per pair.
OUTCOMES per pair (secondary ones in brackets)
  O1 PRIMARY  ensemble mean of the members' days 1-25 mean rainfall, IBERIA
  [O1 for SCANDINAVIA and BRITISH_ISLES; O2 share of members above the control's 90th
   percentile of the days 1-25 mean; O3 mean longest run of wet days (>= 1 mm/day)]
  Forced change D = nudged - control; frequency-only change F = synthetic - control.
STATISTIC  R = sum over pairs of F / sum over pairs of D (share of the forced change
  that regime frequency reproduces); 95% interval from 10,000 resamples of centres
  (pairs of a centre together; synthetic replicates averaged). Also the forced change
  itself with its interval, and the regime-frequency changes.
READING (primary, fixed now)
  If the interval of the forced change D includes 0: "no forced change to explain".
  Else: R interval inside [0.8, 1.2] -> "regime frequency explains the forced rainfall
  change"; R interval entirely below 0.8 -> "rainfall within regimes must change"
  (frequency-only reproduces less than 80%); otherwise "inconclusive".
  Secondary outcomes and k = 6 read the same way, reported without a verdict on the
  primary. Decomposition (descriptive) of O2: sum dpi_k r_k-bar + sum pi_k-bar dr_k.

Output: results/current/8_experiment/snapsi_regime_precip.json
Run:    snapsi_regime_precip.py [--selftest]
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
DATA = HERE.parents[0] / "03_data_ingestion" / "snapsi_regime_precip.parquet"
FIELD = HERE.parents[0] / "03_data_ingestion" / "snapsi_regime_psl.npz"
NAME = "snapsi_regime_precip"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
BOXES = ["IBERIA", "SCANDINAVIA", "BRITISH_ISLES"]
N_EOF, K_PRIMARY, K_SENS, N_RESTART, N_SYN, N_BOOT, WIN = 10, 4, 6, 50, 200, 10000, 5
WET = 1.0


# ------------------------------------------------------------------ regimes
def regimes(con, nud, lat, k, seed):
    """con, nud: (n_days, n_grid) anomaly maps. Returns labels for control and nudged."""
    from sklearn.cluster import KMeans
    w = np.sqrt(np.cos(np.deg2rad(lat)))
    X = con * w
    mu = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    E = Vt[:N_EOF]
    pc_c = (X - mu) @ E.T; pc_n = (nud * w - mu) @ E.T
    km = KMeans(n_clusters=k, n_init=N_RESTART, random_state=seed).fit(pc_c)
    return km.labels_, km.predict(pc_n), {"explained_var_10eof": float((S[:N_EOF] ** 2).sum() / (S ** 2).sum())}


# ------------------------------------------------------------------ construction
def outcomes(rain, members, days):
    """rain: (n_rows, n_box) daily; members/days: row keys. Returns per member days-1-25
    mean and longest wet run, per box."""
    df = pd.DataFrame(rain, columns=BOXES); df["m"] = members; df["d"] = days
    g = df.sort_values(["m", "d"]).groupby("m")
    mean = g[BOXES].mean()
    def run(x):
        best = r = 0
        for v in x >= WET:
            r = r + 1 if v else 0; best = max(best, r)
        return best
    wet = g[BOXES].agg(run)
    return mean, wet


def synth(con_lab, con_day, con_rain, nud_lab, nud_day, rng):
    out = np.empty((len(nud_lab), con_rain.shape[1]))
    for i, (j, d) in enumerate(zip(nud_lab, nud_day)):
        for width in (WIN, 2 * WIN, 99):
            pool = np.flatnonzero((con_lab == j) & (np.abs(con_day - d) <= width))
            if len(pool):
                break
        if not len(pool):
            pool = np.arange(len(con_lab))
        out[i] = con_rain[rng.choice(pool)]
    return out


def pair_stats(p, lab_c, lab_n, rng):
    """p: dict with control/nudged member, day, rain arrays for one pair."""
    c, n = p["control"], p["nudged"]
    mc, wc = outcomes(c["rain"], c["member"], c["day"]); mn, wn = outcomes(n["rain"], n["member"], n["day"])
    q90 = mc.quantile(0.9)
    obs = {"O1": mn.mean() - mc.mean(), "O2": (mn > q90).mean() - (mc > q90).mean(), "O3": wn.mean() - wc.mean()}
    syn = {k: [] for k in obs}
    for _ in range(N_SYN):
        r = synth(lab_c, c["day"], c["rain"], lab_n, n["day"], rng)
        ms, ws = outcomes(r, n["member"], n["day"])
        syn["O1"].append(ms.mean() - mc.mean()); syn["O2"].append((ms > q90).mean() - (mc > q90).mean())
        syn["O3"].append(ws.mean() - wc.mean())
    return {k: (obs[k], pd.concat(v, axis=1).mean(axis=1)) for k, v in syn.items()}


def share(rows, rng):
    """rows: list of (centre, D Series, F Series). R = sum F / sum D per box, centre bootstrap."""
    cs = sorted({r[0] for r in rows})
    def agg(sel):
        D = sum(r[1] for r in sel); F = sum(r[2] for r in sel)
        return D, F, F / D.where(D.abs() > 1e-12)
    D, F, R = agg(rows); bsD, bsR = [], []
    for _ in range(N_BOOT):
        pick = rng.choice(cs, len(cs)); sel = [r for c in pick for r in rows if r[0] == c]
        d_, _, r_ = agg(sel); bsD.append(d_); bsR.append(r_)
    bsD, bsR = pd.concat(bsD, axis=1), pd.concat(bsR, axis=1)
    out = {}
    for b in BOXES:
        n = len(rows)
        dci = np.nanpercentile(bsD.loc[b], [2.5, 97.5]); rci = np.nanpercentile(bsR.loc[b], [2.5, 97.5])
        out[b] = {"forced_change_mean_per_pair": round(float(D[b] / n), 4), "forced_ci95_per_pair": [round(float(x / n), 4) for x in dci],
                  "freq_only_change_mean_per_pair": round(float(F[b] / n), 4), "R": round(float(R[b]), 4),
                  "R_ci95": [round(float(x), 4) for x in rci]}
    return out


def reading(blk):
    lo, hi = blk["forced_ci95_per_pair"]
    if lo <= 0 <= hi:
        return "no forced change to explain"
    a, b = blk["R_ci95"]
    if 0.8 <= a and b <= 1.2:
        return "regime frequency explains the forced rainfall change"
    if b < 0.8:
        return "rainfall within regimes must change"
    return "inconclusive"


def analyse(pairs, fields, rng):
    """pairs: {(centre, init): {arm: {member, day, rain}}}; fields: {centre: {arm: (maps, keys)}}."""
    res = {}
    for k in (K_PRIMARY, K_SENS):
        rows = {o: [] for o in ("O1", "O2", "O3")}; freq = {}
        for centre, F in fields.items():
            lab_c, lab_n, info = regimes(F["control"]["maps"], F["nudged"]["maps"], F["lat"], k, SEED + k)
            kc, kn = F["control"]["keys"], F["nudged"]["keys"]          # (init, member, day) per map row
            lc = {t: l for t, l in zip(map(tuple, kc), lab_c)}; ln = {t: l for t, l in zip(map(tuple, kn), lab_n)}
            freq[centre] = {"control": np.bincount(lab_c, minlength=k).tolist(), "nudged": np.bincount(lab_n, minlength=k).tolist(), **info}
            for (cc, init), p in pairs.items():
                if cc != centre:
                    continue
                labc = np.array([lc.get((init, m, d), -1) for m, d in zip(p["control"]["member"], p["control"]["day"])])
                labn = np.array([ln.get((init, m, d), -1) for m, d in zip(p["nudged"]["member"], p["nudged"]["day"])])
                if (labc < 0).any() or (labn < 0).any():
                    continue                                              # rainfall day without a pressure map
                st = pair_stats(p, labc, labn, rng)
                for o, (D, Fs) in st.items():
                    rows[o].append((centre, D, Fs))
        blk = {o: share(rows[o], rng) for o in rows if rows[o]}
        for o in blk:
            for b in BOXES:
                blk[o][b]["reading"] = reading(blk[o][b])
        blk["regime_counts"] = freq; blk["n_pairs"] = len(rows["O1"])
        res[f"k{k}"] = blk
    res["primary_reading"] = res[f"k{K_PRIMARY}"]["O1"]["IBERIA"]["reading"]
    return res


def selftest():
    """Synthetic: (a) rainfall depends on regime only and forcing changes regime frequency
    -> R ~ 1; (b) forcing adds rain within regimes, frequencies fixed -> R ~ 0."""
    rng = np.random.default_rng(3)
    lat = np.repeat(np.arange(30, 81, 5), 25).astype(float)
    cen = rng.normal(size=(4, lat.size)) * 2
    def make(probs, extra):
        maps, keys, mem, day, rain = [], [], [], [], []
        for m in range(30):
            for d in range(1, 26):
                j = rng.choice(4, p=probs)
                maps.append(cen[j] + rng.normal(size=lat.size)); keys.append(("s1", f"m{m}", d))
                mem.append(f"m{m}"); day.append(d)
                rain.append(np.maximum(0, [1 + j + extra + rng.normal(), 2 - 0.3 * j + rng.normal(), 1.5 + rng.normal()]))
        return np.array(maps), np.array(keys, dtype=object), np.array(mem), np.array(day), np.array(rain)
    out = {}
    for case, (pc, pn, ex) in {"frequency_only": ([0.4, 0.3, 0.2, 0.1], [0.1, 0.2, 0.3, 0.4], 0.0),
                               "within_regime": ([0.25] * 4, [0.25] * 4, 1.0)}.items():
        mc, kc, mmc, dc, rc = make(pc, 0.0); mn, kn, mmn, dn, rn = make(pn, ex)
        pairs = {("X", "s1"): {"control": {"member": mmc, "day": dc, "rain": rc}, "nudged": {"member": mmn, "day": dn, "rain": rn}},
                 ("Y", "s1"): {"control": {"member": mmc, "day": dc, "rain": rc}, "nudged": {"member": mmn, "day": dn, "rain": rn}}}
        fields = {c: {"lat": lat, "control": {"maps": mc, "keys": kc}, "nudged": {"maps": mn, "keys": kn}} for c in ("X", "Y")}
        global N_SYN, N_BOOT
        N_SYN, N_BOOT = 20, 200
        r = analyse(pairs, fields, rng)
        out[case] = {"R_IBERIA": r["k4"]["O1"]["IBERIA"]["R"], "reading": r["k4"]["O1"]["IBERIA"]["reading"]}
    print(json.dumps(out, indent=1))
    ok = out["frequency_only"]["R_IBERIA"] > 0.8 and out["within_regime"]["R_IBERIA"] < 0.3
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    if not DATA.exists() or not FIELD.exists():
        sys.exit("no data: run 03_data_ingestion/acquire_snapsi_regime_precip.py")
    d = pd.read_parquet(DATA); z = np.load(FIELD, allow_pickle=True)
    pairs = {}
    for (c, init, arm), g in d.groupby(["centre", "init", "experiment"]):
        pairs.setdefault((c, init), {})[arm] = {"member": g.member.values, "day": g.day.values, "rain": g[BOXES].values}
    pairs = {k: v for k, v in pairs.items() if len(v) == 2}
    fields = {}
    for c in sorted({k[0] for k in pairs}):
        fields[c] = {"lat": z["lat"], **{arm: {"maps": z[f"{c}|{arm}|maps"], "keys": z[f"{c}|{arm}|keys"]} for arm in ("control", "nudged")}}
    res = {"registered": "2026-10-08, committed before any value was fetched", "seed": SEED, "n_syn": N_SYN, "n_boot": N_BOOT,
           **analyse(pairs, fields, np.random.default_rng(SEED))}
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(json.dumps({k: v for k, v in res.items() if k != "k6"}, indent=1, default=float)[:4000])
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
