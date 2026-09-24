#!/usr/bin/env python3
"""
within_model_check.py
=====================
IS THE HEADLINE R^2 WITHIN-MODEL SKILL, OR BETWEEN-MODEL SIGNAL?

THE PROBLEM, FOUND BY STRESS-TESTING RATHER THAN BY A REFEREE
  `headline_stress_test.py` computed CV R^2 separately inside each ensemble
  member and found it scattered around ZERO (-0.18 to +0.17 across 20 members),
  while the pooled cross-member estimate is +0.115.

  That gap has an innocent reading and a fatal one.

  INNOCENT: within-member n is only ~94 against 38 features, so per-member CV is
  swamped by estimation noise and its near-zero average means nothing.

  FATAL: GroupKFold by member stops the model memorising a member, but it does
  NOT stop it exploiting BETWEEN-MODEL structure. If models differ systematically
  in both their typical vortex anomaly and their typical surface response, a
  regression can learn that cross-model relationship and score well on held-out
  members without any ability to rank events INSIDE one climate. The field's
  question -- given an SSW in the real atmosphere, will this one couple
  downward? -- is a within-climate question. Between-model skill does not answer
  it, and reporting 0.115 as if it did would be wrong.

THE TEST
  Standardise Y and every feature WITHIN each member (subtract that member's mean,
  divide by its sd) before pooling. That removes all between-model level
  differences by construction while leaving every within-model event-to-event
  relationship intact. Then rerun the identical pooled GroupKFold pipeline.

    R^2 holds near 0.115  ->  the skill is WITHIN-model. Headline stands.
    R^2 collapses to ~0   ->  the skill was BETWEEN-model. The headline must be
                              restated as a bound and the paper's central number
                              changes.

  Pseudo-events are carried through the same transformation, so the artefact
  share measured in Attack 1 can be re-read on the same footing.

SIZE-MATCHED NULL AND INTERVALS (added 2026-09-23)
  The original pseudo arm pools 6 draws per member: 11,321 pseudo-events
  (11,323 in the 2026-08 environment) against 1,888 real ones. With 38 features a ridge fit on 6x the data scores a
  higher out-of-sample R^2 for that reason alone, which inflates the pseudo
  baseline and biases "event-specific = within - pseudo_within" toward zero.
  So the comparison is repeated against K_DRAWS independent pseudo draws, EACH
  the same size as the real set (same members, same calendar days, drawn from
  event-cleaned days -- the existing `draw_pseudo`). Each draw's within-model CV
  R^2 is one sample from the no-SSW distribution at the real sample size:
    event_specific_matched = within_real - mean(null)
    p = (1 + #{null >= within_real}) / (K + 1)
  A second, independent interval comes from a member-cluster bootstrap
  (N_BOOT resamples of the 20 members with replacement; duplicates keep one
  group label so they are held out together, and one pseudo draw per member
  per replicate so duplicates are identical in both arms).
  The original 6-draw numbers are still computed first, as the reproduction
  check, with the original seed and rng sequence untouched.

PARALLEL, AND STILL DETERMINISTIC (2026-09-24)
  The matched draws (per member), the null fits (per draw) and the bootstrap
  replicates run on N_WORKERS processes. Every task draws from its OWN stream,
  np.random.SeedSequence(seed).spawn(n)[i], fixed by the task's index -- never
  from a shared generator -- so the result does not depend on how tasks are
  scheduled or how many workers there are. The bootstrap uses ONE resample per
  replicate for all three tiers (paired across tiers). BLAS is pinned to one
  thread per worker so the processes do not oversubscribe the cores.

Output: within_model_check.json
"""
import json
import os
import sys
import warnings
import zlib

# One BLAS thread per process, set before numpy loads: the parallel sections
# below run N_WORKERS processes and would otherwise oversubscribe every core.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "6_predictability"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402
import predictability_ceiling as P                  # noqa: E402
import headline_stress_test as H                    # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
SEED = 20260805
K_DRAWS = 1000        # size-matched pseudo draws per member
N_BOOT = 1000         # member-cluster bootstrap replicates
SEED_MATCHED = zlib.crc32(b"within_model_check_size_matched") % (2 ** 32)
SEED_BOOT = zlib.crc32(b"within_model_check_cluster_bootstrap") % (2 ** 32)
TIERS = ((1, "P1 pre-onset"), (2, "P2 at-onset"), (3, "P3 + post-onset stratosphere"))
N_WORKERS = max(1, min(20, (os.cpu_count() or 2) - 2))


def demean_within(X, y, g):
    """Remove each member's own mean and scale, leaving only within-member variation."""
    Xo = X.copy().astype(float)
    yo = y.copy().astype(float)
    for m in np.unique(g):
        i = g == m
        for c in range(Xo.shape[1]):
            v = Xo[i, c]
            s = np.nanstd(v)
            Xo[i, c] = (v - np.nanmean(v)) / (s if s > 1e-12 else 1.0)
        sy = np.nanstd(yo[i])
        yo[i] = (yo[i] - np.nanmean(yo[i])) / (sy if sy > 1e-12 else 1.0)
    return Xo, yo


def prepare_member(f):
    """Everything one member contributes, or None if it is not usable.

    Shared by the sequential real/6-draw build and the parallel matched draws,
    so both see exactly the same onsets, clean days and climatology.
    """
    try:
        m = EP.load_member(f)
    except Exception:
        return None
    m2 = m[np.isin(m.index.month, SEASON)].dropna()
    if len(m2) < 2000:
        return None
    on = EP.detect_ssw(m2["u10"].values, m2.index)
    if len(on) < 15:
        return None
    am = m2["am"]
    msk = C6.influence_mask(am.index, on)
    cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
    cln = am[~msk].index
    dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
    ds = xr.open_dataset(f)
    lat, plev, u = ds["lat"].values, ds["plev"].values, ds["u_zm"].values
    ds.close()
    idx = m.index
    if len(idx) != u.shape[0]:
        return None
    return dict(on=on, am=am, cl=cl, cln=cln, dy=dy, lat=lat, plev=plev, u=u, idx=idx)


def matched_draws(task):
    """K size-matched pseudo draws for one member, from that member's own stream."""
    i, f = task
    c = prepare_member(Path(f))
    rng = np.random.default_rng(np.random.SeedSequence(SEED_MATCHED).spawn(i + 1)[i])
    n_real = len(C6.anom(c["am"], c["on"], c["cl"], OUT_WIN))
    Fs, Ys, cols = [], [], None
    for _ in range(K_DRAWS):
        p_on = C6.draw_pseudo(c["cln"], c["dy"], rng)
        Fk = H.features_for(c["u"], c["plev"], c["lat"], c["idx"], p_on)
        Yk = C6.anom(c["am"], p_on, c["cl"], OUT_WIN)
        if len(Yk) < n_real:
            pad = n_real - len(Yk)
            Fk = pd.concat([Fk, pd.DataFrame(np.nan, index=range(pad),
                                             columns=Fk.columns)], ignore_index=True)
            Yk = np.concatenate([Yk, np.full(pad, np.nan)])
        Fs.append(Fk.values); Ys.append(Yk); cols = list(Fk.columns)
    return Path(f).stem, np.stack(Fs), np.stack(Ys), cols


# Shared read-only state for the forked workers of the null and the bootstrap.
_S = {}


def _stack(ks, ci, mem_list):
    M = _S["matched"]
    X = np.concatenate([M[m][0][ks[m]][:, ci] for m in mem_list])
    y = np.concatenate([M[m][1][ks[m]] for m in mem_list])
    g = np.concatenate([np.full(M[m][1].shape[1], m) for m in mem_list])
    return X, y, g


def null_task(k):
    """Within-model CV R^2 of matched draw k, all members, every tier."""
    out = []
    for ci in _S["tier_ci"]:
        X, y, g = _stack({m: k for m in _S["members"]}, ci, _S["members"])
        Xn, yn = demean_within(X, y, g)
        out.append(P.cv_r2(Xn, yn, g, "ridge"))
    return out


def boot_task(b):
    """One member-cluster resample, used for all tiers: real minus pseudo."""
    members = _S["members"]
    rng = np.random.default_rng(np.random.SeedSequence(SEED_BOOT).spawn(b + 1)[b])
    pick = list(rng.choice(members, size=len(members), replace=True))
    # sorted(): iterating a set of strings follows the per-process string hash.
    ks = {m: int(rng.integers(K_DRAWS)) for m in sorted(set(pick))}
    ir = np.concatenate([_S["real_by"][m] for m in pick])
    gb = np.concatenate([np.full(len(_S["real_by"][m]), m) for m in pick])
    out = []
    for ci, cols in zip(_S["tier_ci"], _S["tier_cols"]):
        Xb, yb = demean_within(_S["Xr"][cols].values[ir], _S["yr"][ir], gb)
        rb = P.cv_r2(Xb, yb, gb, "ridge")
        Xp_, yp_, gp_ = _stack(ks, ci, pick)
        Xpn, ypn = demean_within(Xp_, yp_, gp_)
        pb = P.cv_r2(Xpn, ypn, gp_, "ridge")
        out.append(rb - pb if np.isfinite(rb) and np.isfinite(pb) else np.nan)
    return out


def main():
    rng = np.random.default_rng(SEED)
    real_rows, real_Y, real_g = [], [], []
    ps_rows, ps_Y, ps_g = [], [], []
    used = []

    print("rebuilding CMIP6 real and pseudo tables ...", flush=True)
    for f in sorted(RAW.glob("*_zm.nc")):
        c = prepare_member(f)
        if c is None:
            continue
        used.append(f)
        on, am, cl, cln, dy = c["on"], c["am"], c["cl"], c["cln"], c["dy"]
        lat, plev, u, idx = c["lat"], c["plev"], c["u"], c["idx"]
        Yr = C6.anom(am, on, cl, OUT_WIN)
        real_rows.append(H.features_for(u, plev, lat, idx, on))
        real_Y.append(Yr)
        real_g.append(np.full(len(Yr), f.stem))
        for _ in range(6):
            p_on = C6.draw_pseudo(cln, dy, rng)
            if len(p_on) < 10:
                continue
            Yp = C6.anom(am, p_on, cl, OUT_WIN)
            ps_rows.append(H.features_for(u, plev, lat, idx, p_on))
            ps_Y.append(Yp)
            ps_g.append(np.full(len(Yp), f.stem))
        del c, u

    # K size-matched draws per member, each with the real set's onset count. A
    # draw can still carry FEWER USABLE rows than the real set: a pseudo onset
    # whose +8..+52 window lacks enough in-season days gets a NaN outcome
    # (C6.anom), which cv_r2 drops. Those are counted below as "short". The NaN
    # padding is a guard in case draw_pseudo ever returns fewer onsets (a
    # day-of-year with no clean day); it did not fire on 2026-09-23.
    ctx = mp.get_context("fork")
    print(f"generating {K_DRAWS} matched draws x {len(used)} members "
          f"on {N_WORKERS} processes ...", flush=True)
    with ctx.Pool(N_WORKERS) as pool:
        matched = {stem: (F, Y, cols) for stem, F, Y, cols in
                   pool.map(matched_draws, [(i, str(f)) for i, f in enumerate(used)])}

    Xr = pd.concat(real_rows, ignore_index=True)
    yr = np.concatenate(real_Y)
    gr = np.concatenate(real_g)
    Xp = pd.concat(ps_rows, ignore_index=True)
    yp = np.concatenate(ps_Y)
    gp = np.concatenate(ps_g)
    print(f"real {np.isfinite(yr).sum():,} | pseudo {np.isfinite(yp).sum():,} | "
          f"{len(np.unique(gr))} members")

    res = {}
    print("\n" + "=" * 82)
    print("=== POOLED (as published) vs WITHIN-MODEL STANDARDISED ===")
    print(f"{'tier':<32s} {'pooled':>9s} {'within':>9s} {'retained':>9s} "
          f"{'pseudo-w':>9s}")
    print("-" * 76)
    for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                      (3, "P3 + post-onset stratosphere")):
        cols = P.tier_cols(Xr.columns, tier) + ["doy_sin", "doy_cos"]
        Xv, Xpv = Xr[cols].values, Xp[cols].values
        pooled = P.cv_r2(Xv, yr, gr, "ridge")
        Xw, yw = demean_within(Xv, yr, gr)
        within = P.cv_r2(Xw, yw, gr, "ridge")
        Xpw, ypw = demean_within(Xpv, yp, gp)
        pseudo_w = P.cv_r2(Xpw, ypw, gp, "ridge")
        ret = within / pooled if pooled else np.nan
        res[lab] = {"pooled_cv_r2": round(pooled, 4),
                    "within_model_cv_r2": round(within, 4),
                    "retained_fraction": None if not np.isfinite(ret) else round(float(ret), 3),
                    "pseudo_within_cv_r2": round(pseudo_w, 4),
                    "event_specific": round(float(within - pseudo_w), 4)}
        print(f"{lab:<32s} {pooled:+9.4f} {within:+9.4f} {ret:9.2f} {pseudo_w:+9.4f}")

    # The verdict this block used to print ("Headline stands") tested
    # within > 0.3 x pooled, which is not the question. The quantity that
    # matters is within minus pseudo-within, reported below with intervals.
    print(f"\n  event-specific within-model share, 6-draw pooled pseudo arm (original):")
    for k, v in res.items():
        print(f"    {k:<32s} {v['event_specific']:+.4f}")

    # ------------------------------------------- size-matched null + bootstrap
    members = sorted(matched)
    real_by = {m: np.flatnonzero(gr == m) for m in members}
    cols_all = matched[members[0]][2]
    assert list(Xr.columns) == cols_all, "real and matched feature columns differ"
    sizes = {m: [int(np.isfinite(matched[m][1][k]).sum()) for k in range(K_DRAWS)]
             for m in members}
    # Real outcomes are all finite (1,888 of 1,888), so len(real_by[m]) is the
    # real arm's usable count; a draw is short when fewer of its outcomes are.
    short = sum(int(np.sum(np.array(v) < len(real_by[m]))) for m, v in sizes.items())
    print(f"\n=== SIZE-MATCHED NULL: {K_DRAWS} draws, each the size of the real set ===")
    print(f"  draws with fewer usable outcomes than the real set (NaN outcome "
          f"window): "
          f"{short} of {K_DRAWS * len(members)} member-draws")

    tier_cols = [P.tier_cols(Xr.columns, t) + ["doy_sin", "doy_cos"] for t, _ in TIERS]
    _S.update(matched=matched, members=members, real_by=real_by, Xr=Xr, yr=yr,
              tier_cols=tier_cols,
              tier_ci=[[cols_all.index(c) for c in cols] for cols in tier_cols])
    print(f"null: {K_DRAWS} draws; bootstrap: {N_BOOT} replicates; "
          f"{N_WORKERS} processes ...", flush=True)
    with ctx.Pool(N_WORKERS) as pool:
        nulls = np.array(pool.map(null_task, range(K_DRAWS), chunksize=10), float)
        boots = np.array(pool.map(boot_task, range(N_BOOT), chunksize=10), float)

    for ti, (tier, lab) in enumerate(TIERS):
        cols = tier_cols[ti]
        Xw, yw = demean_within(Xr[cols].values, yr, gr)
        real_w = P.cv_r2(Xw, yw, gr, "ridge")
        null = nulls[:, ti]
        null = null[np.isfinite(null)]
        es = real_w - float(null.mean())
        p_null = (1 + int(np.sum(null >= real_w))) / (len(null) + 1)
        diffs = boots[:, ti]
        diffs = diffs[np.isfinite(diffs)]
        res[lab]["size_matched"] = {
            "within_real_cv_r2": round(float(real_w), 4),
            "null_mean": round(float(null.mean()), 4),
            "null_sd": round(float(null.std(ddof=1)), 4),
            "null_q025_q975": [round(float(np.quantile(null, 0.025)), 4),
                               round(float(np.quantile(null, 0.975)), 4)],
            "event_specific": round(es, 4),
            "event_specific_CI95_from_null": [
                round(float(real_w - np.quantile(null, 0.975)), 4),
                round(float(real_w - np.quantile(null, 0.025)), 4)],
            "p_null_ge_real": round(p_null, 4),
            "n_null_draws": int(len(null)),
            "bootstrap_diff_mean": round(float(diffs.mean()), 4),
            "bootstrap_CI95": [round(float(np.quantile(diffs, 0.025)), 4),
                               round(float(np.quantile(diffs, 0.975)), 4)],
            "n_boot_valid": int(len(diffs)),
            "size_bias_of_6draw_pooled_arm": round(
                res[lab]["pseudo_within_cv_r2"] - float(null.mean()), 4)}
        sm = res[lab]["size_matched"]
        print(f"  {lab:<30s} real {real_w:+.4f} | null {sm['null_mean']:+.4f} "
              f"(sd {sm['null_sd']:.4f}) | event-specific {es:+.4f} "
              f"[{sm['event_specific_CI95_from_null'][0]:+.4f}, "
              f"{sm['event_specific_CI95_from_null'][1]:+.4f}] p={p_null:.4f} | "
              f"boot [{sm['bootstrap_CI95'][0]:+.4f}, {sm['bootstrap_CI95'][1]:+.4f}] | "
              f"6-draw size bias {sm['size_bias_of_6draw_pooled_arm']:+.4f}", flush=True)

    import sklearn
    res["_meta"] = {
        "seed_original": SEED, "seed_size_matched": SEED_MATCHED,
        "seed_bootstrap": SEED_BOOT, "n_workers": N_WORKERS,
        "streams": "SeedSequence(seed).spawn(n)[task index]; scheduling-independent",
        "k_draws": K_DRAWS, "n_boot": N_BOOT,
        "n_real_events": int(np.isfinite(yr).sum()),
        "n_pseudo_events_6draw": int(np.isfinite(yp).sum()),
        "n_members": len(members), "outcome_window_days": list(OUT_WIN),
        "season_months": list(SEASON), "input_dir": str(RAW.relative_to(ROOT)),
        "matched_draws_short": short,
        "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                     "sklearn": sklearn.__version__, "xarray": xr.__version__}}
    (RESULTS / "within_model_check.json").write_text(json.dumps(res, indent=2),
                                                  encoding="utf8", newline="\n")
    print("\nSaved -> within_model_check.json")


if __name__ == "__main__":
    main()
