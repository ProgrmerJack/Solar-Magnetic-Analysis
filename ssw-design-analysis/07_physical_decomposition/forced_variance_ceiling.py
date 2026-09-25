#!/usr/bin/env python3
"""
forced_variance_ceiling.py
==========================
THE POSITIVE CLAIM UNDERNEATH THE WHOLE PROJECT, AND A HARD CEILING ON EVERY
CLASSIFIER IN THIS LITERATURE.

THE QUESTION NOBODY HAS ASKED
  The field asks "which SSWs propagate downward?" and builds classifiers to find
  out. That question presupposes an answer to a prior one: DO SSWs DIFFER FROM ONE
  ANOTHER IN THEIR FORCED SURFACE RESPONSE AT ALL? If every SSW forces the same
  expected surface anomaly and the observed spread between events is internal
  tropospheric noise, then there is nothing to classify, and every classifier --
  however clever, however physical -- is sorting noise.

  `is_downward_propagation_a_class.py` established that the per-event response
  distribution is a PURE TRANSLATION of the ordinary winter distribution. That
  already implies this, because a translation changes location and not spread. But
  it was never stated as a variance claim and never given a number. This does both.

THE DECOMPOSITION
  Write the per-event surface anomaly as

      Y_i  =  mu + S + f_i + eps_i           at real SSWs
      Y_j  =  mu     +       eps_j           at event-free pseudo-onsets

  S      the mean forced shift, common to every event -- the real, large SSW effect
  f_i    THE EVENT-SPECIFIC part of the forced response, E[f] = 0, Var(f) = s_f^2
  eps    internal tropospheric noise, Var(eps) = s_e^2

  With f independent of eps,

      Var(Y | SSW)     =  s_f^2 + s_e^2
      Var(Y | pseudo)  =          s_e^2
      ------------------------------------------------
      s_f^2  =  Var(Y | SSW)  -  Var(Y | pseudo)                          (*)

  s_f is the ONLY quantity a classifier can recover. Everything else is noise that
  no predictor, in principle, can see.

THE CEILING, WHICH IS THE POINT -- STATED CAREFULLY
  Suppose a PERFECT classifier: one that knows f_i exactly for every event. Split
  the events at the q-quantile of f and difference the group means. For Gaussian f

      contrast(q)  =  s_f * phi(z_q) * (1/q + 1/(1-q))                    (**)

  A FIRST DRAFT OF THIS SCRIPT USED 2*sqrt(2/pi) = 1.596, the median-split value,
  and called it the maximum. THAT IS WRONG, AND IN THE DIRECTION THAT FLATTERS THE
  ARGUMENT. (**) is MINIMISED at q = 0.5 and grows without bound as the split gets
  more extreme -- 1.596 at q=0.5, 1.655 at q=0.7, 1.950 at q=0.9. A classifier free
  to choose its split point has no finite ceiling at all, because it can always
  compare the single most extreme event against the rest.

  So the honest quantity is the maximum contrast achievable AT THE SPLIT FRACTION
  THE CRITERION ITSELF PRODUCES. Each contrast is therefore computed here with the
  published criterion, on the same outcome as s_f, and (**) is evaluated at the
  DW fraction it gives. For the split
  fractions actually in use (54-70% DW) the factor lands in 1.60-1.66, so the
  numbers barely move -- but the claim is now the one the algebra supports.

ASSUMPTION, STATED
  (*) needs f independent of eps. If a large forced response also suppressed
  internal variance, the decomposition would misattribute. Nothing here tests
  that, and it is recorded as an assumption rather than buried.

WHY THE TEST IS CONSERVATIVE AGAINST ITS OWN CONCLUSION
  Pseudo-onsets are drawn only from days more than 135 days from every real
  onset, so no pseudo-event's -60..+75 day neighbourhood overlaps a real event's.
  That excludes the most disturbed periods of winter, which DEFLATES
  Var(Y | pseudo) and therefore INFLATES s_f^2 by (*). The bias runs against the
  finding of small s_f, not toward it.

THREE CONTROLS, BECAUSE (*) HAS ASSUMPTIONS
  C1 PRE-ONSET PLACEBO. Run the identical decomposition on the window BEFORE onset
     (days -52..-8), where by construction s_f = 0. A non-zero estimate there is
     the method's own bias and is subtracted from nothing -- it is reported as the
     error bar the estimate has to beat.
  C2 KNOWN TRUTH. Simulate with s_f injected at known values and check recovery.
     An estimator this project has not validated under known truth is not used.
  C3 NOISE-EQUALITY. (*) assumes s_e^2 is the same at events and pseudo-events. If
     an SSW changes tropospheric variance itself, that assumption fails. Tested by
     comparing pre-onset variance between the two sets.

Output: forced_variance_ceiling.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "6_predictability"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "08_literature_audit"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import selection_on_outcome as SO                   # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402
import predictability_ceiling as PC                 # noqa: E402
import recompute_published_criterion as RPC         # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
OUT_WIN = (8, 52)
PRE_WIN = (-52, -8)          # the placebo window, same length, before onset
SEASON = (11, 12, 1, 2, 3, 4)
N_PSEUDO = 400               # pseudo-onset draws (observations)
N_BOOT = 3000
SEED = 20260803

# The DW-minus-NDW contrast each ceiling is compared with is COMPUTED here, on the
# same per-event outcome the ceiling's s_f comes from, with the split fraction q
# that the criterion produced on those events. Until 2026-09-25 this was a typed-in
# list: ERA5 NAM contrasts were tagged as matched to the CPC-AO s_f, and the CMIP6
# row was a median split (q = 0.5 by construction) from an older event set. The
# cross-index rows are gone; their matched implied R^2 is in their own scripts.


def matched_contrast(Y, lab):
    """(contrast DW - NDW, DW fraction, n) on events with a finite outcome and label."""
    Y = np.asarray(Y, float)
    lab = np.asarray(lab, float)
    ok = np.isfinite(Y) & np.isfinite(lab)
    g = lab[ok].astype(bool)
    return float(Y[ok][g].mean() - Y[ok][~g].mean()), float(g.mean()), int(ok.sum())


def sigma_f2(Y_ev, Y_ps):
    """The estimator (*): excess event variance over matched pseudo-event variance."""
    a = Y_ev[np.isfinite(Y_ev)]
    b = Y_ps[np.isfinite(Y_ps)]
    return float(np.var(a, ddof=1) - np.var(b, ddof=1))


def split_factor(q):
    """(**) contrast per unit s_f for a perfect classifier splitting at quantile q.

    phi(z_q) * (1/q + 1/(1-q)). Minimised at q=0.5 (=1.596), rising for any more
    extreme split -- which is why the split fraction must come from the study.
    """
    from scipy.stats import norm
    q = min(max(float(q), 1e-3), 1 - 1e-3)
    return float(norm.pdf(norm.ppf(q)) * (1.0 / q + 1.0 / (1.0 - q)))


def ceiling(s2, q=0.5):
    """Max DW-NDW contrast a PERFECT classifier could produce at split fraction q."""
    return split_factor(q) * np.sqrt(s2) if s2 > 0 else 0.0


def winter_of(idx):
    t = pd.DatetimeIndex(idx)
    return np.where(t.month >= 11, t.year + 1, t.year)


def boot_paired(Y_post, Y_pre, wid, pool_post, pool_pre, rng, n_boot=N_BOOT):
    """PAIRED bootstrap of the placebo-corrected estimate.

    The placebo window is not decoration. In CMIP6 the pre-onset window -- where
    the event cannot possibly have created event-to-event differences -- returns
    s_f^2 = +0.067 against +0.075 post-onset. Nearly the whole apparent forced
    variance is there BEFORE the event. The cause is that SSWs are not randomly
    timed: they cluster in disturbed winters whose tropospheric variance is
    already elevated (the C3 ratio is 1.14), and pseudo-onsets matched only on
    day-of-year do not reproduce that. Reporting the uncorrected number as
    "event-specific forced response" would have been a straight false positive.

    Post and pre are computed on THE SAME resampled winters so the difference is
    paired and the shared selection confound cancels.

    HOW TO READ THE PAIR, which matters more than either number alone:
      uncorrected  -> UPPER bound on s_f. It still contains the confound.
      corrected    -> LOWER bound on s_f. If preconditioning genuinely makes some
                      events force a larger response, part of the pre-onset excess
                      is signal and subtracting it removes real effect.
    The ceiling is therefore quoted from the UNCORRECTED upper bound, which is the
    choice that works against the conclusion being argued for.
    """
    uw = np.unique(wid)
    fin_post = pool_post[np.isfinite(pool_post)]
    fin_pre = pool_pre[np.isfinite(pool_pre)]
    out = []
    for _ in range(n_boot):
        pick = rng.choice(uw, len(uw), replace=True)
        sel = np.concatenate([np.flatnonzero(wid == q) for q in pick])
        a, b = Y_post[sel], Y_pre[sel]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) < 8 or len(b) < 8:
            continue
        pa = rng.choice(fin_post, len(fin_post), replace=True)
        pb = rng.choice(fin_pre, len(fin_pre), replace=True)
        out.append((np.var(a, ddof=1) - np.var(pa, ddof=1))
                   - (np.var(b, ddof=1) - np.var(pb, ddof=1)))
    return np.array(out)


def boot_ci(Y_ev, wid, pool, rng, n_boot=N_BOOT):
    """Winter-block bootstrap of s_f^2 and of the implied ceiling.

    Events are resampled by WINTER, not individually: two SSWs in the same winter
    share the same seasonal state and are not independent draws.
    """
    uw = np.unique(wid)
    vp = np.var(pool[np.isfinite(pool)], ddof=1)
    out = []
    for _ in range(n_boot):
        pick = rng.choice(uw, len(uw), replace=True)
        sel = np.concatenate([np.flatnonzero(wid == q) for q in pick])
        y = Y_ev[sel]
        y = y[np.isfinite(y)]
        if len(y) < 8:
            continue
        # resample the pseudo pool too, so its sampling error is carried
        pb = rng.choice(pool[np.isfinite(pool)], len(pool[np.isfinite(pool)]),
                        replace=True)
        out.append(np.var(y, ddof=1) - np.var(pb, ddof=1))
    out = np.array(out)
    return out, vp


def validate_known_truth(pool, n_ev, rng, n_rep=400):
    """C2. Inject a KNOWN s_f and check the estimator returns it."""
    print("\n=== C2  RECOVERY UNDER KNOWN TRUTH ===")
    print(f"  n_events = {n_ev}, noise drawn from the pseudo-event pool")
    print(f"  {'true s_f':>9s} {'recovered s_f':>15s} {'bias':>8s} {'sd':>8s}")
    print("  " + "-" * 44)
    base = pool[np.isfinite(pool)]
    out = {}
    for s_true in (0.0, 0.2, 0.4, 0.6, 0.8):
        est = []
        for _ in range(n_rep):
            eps = rng.choice(base, n_ev, replace=True)
            f = rng.normal(0, s_true, n_ev)
            ps = rng.choice(base, len(base), replace=True)
            v = np.var(eps + f, ddof=1) - np.var(ps, ddof=1)
            est.append(np.sqrt(v) if v > 0 else 0.0)
        est = np.array(est)
        out[f"s_f={s_true}"] = {"recovered_mean": round(float(est.mean()), 4),
                               "bias": round(float(est.mean() - s_true), 4),
                               "sd": round(float(est.std()), 4)}
        print(f"  {s_true:9.2f} {est.mean():15.3f} {est.mean() - s_true:+8.3f} "
              f"{est.std():8.3f}")
    return out


def analyse(name, ao_like, onsets, clim, clean_idx, doys, per_event, rng,
            draw, n_pseudo=N_PSEUDO, reported=()):
    print("\n" + "=" * 74)
    print(f"=== {name} ===")
    Y_ev = per_event(ao_like, onsets, clim, OUT_WIN)
    Y_pre = per_event(ao_like, onsets, clim, PRE_WIN)

    ps_out, ps_pre = [], []
    for _ in range(n_pseudo):
        f = draw(clean_idx, doys, rng)
        if len(f) < 8:
            continue
        ps_out.append(per_event(ao_like, f, clim, OUT_WIN))
        ps_pre.append(per_event(ao_like, f, clim, PRE_WIN))
    pool_out = np.concatenate(ps_out)
    pool_pre = np.concatenate(ps_pre)

    ne = int(np.isfinite(Y_ev).sum())
    print(f"  events {ne}, pseudo-events {int(np.isfinite(pool_out).sum()):,}")
    print(f"  mean shift S           = {np.nanmean(Y_ev) - np.nanmean(pool_out):+.3f} sigma")
    print(f"  Var(Y | SSW)           = {np.nanvar(Y_ev, ddof=1):.4f}")
    print(f"  Var(Y | pseudo)        = {np.nanvar(pool_out, ddof=1):.4f}")

    s2 = sigma_f2(Y_ev, pool_out)
    wid = winter_of(onsets)
    bs, _ = boot_ci(Y_ev, wid, pool_out, rng)
    ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    s_pt = np.sqrt(s2) if s2 > 0 else 0.0
    s_hi = np.sqrt(ci[1]) if ci[1] > 0 else 0.0

    print(f"\n  s_f^2 = Var(SSW) - Var(pseudo) = {s2:+.4f}  "
          f"95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]")
    print(f"  s_f   = {s_pt:.3f} sigma   (upper 95% bound {s_hi:.3f} sigma)")

    # C1 placebo
    s2_pre = sigma_f2(Y_pre, pool_pre)
    bs_pre, _ = boot_ci(Y_pre, wid, pool_pre, rng)
    ci_pre = [float(np.percentile(bs_pre, 2.5)), float(np.percentile(bs_pre, 97.5))]
    print(f"\n  C1 PLACEBO, pre-onset days {PRE_WIN[0]}..{PRE_WIN[1]} "
          f"(true s_f = 0 by construction):")
    print(f"     s_f^2 = {s2_pre:+.4f}  95% CI [{ci_pre[0]:+.4f}, {ci_pre[1]:+.4f}]"
          f"  -> {'consistent with zero' if ci_pre[0] <= 0 <= ci_pre[1] else 'NON-ZERO: a selection confound is present and is corrected below'}")

    # placebo-corrected, paired
    bsd = boot_paired(Y_ev, Y_pre, wid, pool_out, pool_pre, rng)
    s2_c = s2 - s2_pre
    ci_c = [float(np.percentile(bsd, 2.5)), float(np.percentile(bsd, 97.5))]
    print(f"     placebo-CORRECTED s_f^2 = {s2_c:+.4f} "
          f"[{ci_c[0]:+.4f}, {ci_c[1]:+.4f}]  "
          f"-> s_f = {np.sqrt(max(s2_c, 0)):.3f} "
          f"(upper 95% {np.sqrt(max(ci_c[1], 0)):.3f})")

    # C3 noise equality
    vr = float(np.nanvar(Y_pre, ddof=1) / np.nanvar(pool_pre, ddof=1))
    print(f"  C3 NOISE EQUALITY, pre-onset variance ratio events/pseudo = {vr:.3f}"
          f"  -> {'ok' if 0.7 <= vr <= 1.4 else 'ASSUMPTION SUSPECT'}")

    # the ceiling, evaluated at each study's OWN split fraction
    print(f"\n  === CEILING ON ANY CLASSIFIER (uncorrected s_f, the conservative arm) ===")
    print(f"  {'reported contrast':<50s} {'q':>5s} {'|obs|':>7s} {'ceil95':>7s} "
          f"{'ratio':>7s}  match")
    print("  " + "-" * 84)
    comp = {}
    for lab, v, q, n in reported:
        v = abs(v)
        c_hi_q = ceiling(ci[1], q)
        ratio = v / c_hi_q if c_hi_q > 0 else np.inf
        comp[lab] = {"contrast_abs": round(v, 4), "split_fraction": round(q, 3),
                     "n_events": n,
                     "ceiling_upper95_at_q": round(float(c_hi_q), 4),
                     "over_upper_bound_ratio":
                         None if not np.isfinite(ratio) else round(float(ratio), 2),
                     "matched_system": True}
        rs = "inf" if not np.isfinite(ratio) else f"{ratio:.1f}x"
        print(f"  {lab:<50s} {q:5.2f} {v:7.3f} {c_hi_q:7.3f} {rs:>7s}  YES")
    c_pt, c_hi = ceiling(s2, 0.5), ceiling(ci[1], 0.5)
    print(f"\n  (median-split reference: point {c_pt:.3f}, upper 95% {c_hi:.3f} sigma)")

    return {"n_events": ne, "n_pseudo": int(np.isfinite(pool_out).sum()),
            "mean_shift": round(float(np.nanmean(Y_ev) - np.nanmean(pool_out)), 4),
            "var_event": round(float(np.nanvar(Y_ev, ddof=1)), 4),
            "var_pseudo": round(float(np.nanvar(pool_out, ddof=1)), 4),
            "sigma_f2": round(s2, 5), "sigma_f2_CI95": [round(c, 5) for c in ci],
            "sigma_f": round(float(s_pt), 4),
            "sigma_f_upper95": round(float(s_hi), 4),
            "placebo_sigma_f2": round(s2_pre, 5),
            "placebo_CI95": [round(c, 5) for c in ci_pre],
            "placebo_clean": bool(ci_pre[0] <= 0 <= ci_pre[1]),
            "corrected_sigma_f2": round(float(s2_c), 5),
            "corrected_CI95": [round(c, 5) for c in ci_c],
            "corrected_sigma_f": round(float(np.sqrt(max(s2_c, 0))), 4),
            "corrected_sigma_f_upper95": round(float(np.sqrt(max(ci_c[1], 0))), 4),
            "pre_onset_variance_ratio": round(vr, 4),
            "noise_equality_ok": bool(0.7 <= vr <= 1.4),
            "ceiling_point": round(float(c_pt), 4),
            "ceiling_upper95": round(float(c_hi), 4),
            "reported_vs_ceiling": comp}, pool_out


def main():
    rng = np.random.default_rng(SEED)
    res = {"outcome_window": list(OUT_WIN), "placebo_window": list(PRE_WIN),
           "results": {}}

    # ------------------------------------------------ observations
    ao = M.load("ao")["y"]
    real = load_catalogue("primary")
    real = real[(real >= ao.index.min()) & (real <= ao.index.max())]
    mask = G.real_influence_mask(ao.index, real)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    clean_idx = G.zone_free_index(ao.index, real)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])

    def _win_anom(s, on, cl, win):
        """Windowed day-of-year anomaly. Must agree with SO.per_event at OUT_WIN --
        asserted below rather than assumed, since the placebo control is only
        meaningful if both windows use an identical estimator."""
        out = []
        for o in pd.DatetimeIndex(on):
            lo = o + pd.Timedelta(days=win[0])
            hi = o + pd.Timedelta(days=win[1])
            v = s[(s.index >= lo) & (s.index <= hi)]
            if len(v) < (win[1] - win[0]) // 2:
                out.append(np.nan)
                continue
            c = cl.reindex(v.index.dayofyear).values
            out.append(float(np.nanmean(v.values - c)))
        return np.array(out)

    ref = SO.per_event(ao, real, clim)
    mine = _win_anom(ao, real, clim, OUT_WIN)
    both = np.isfinite(ref) & np.isfinite(mine)
    dmax = float(np.abs(ref[both] - mine[both]).max())
    print(f"estimator agreement with selection_on_outcome.per_event: "
          f"max|diff| = {dmax:.2e} on {both.sum()} events")
    assert dmax < 1e-9, "windowed estimator disagrees with the project's own"

    # Karpechko et al. (2017) conditions 1-3 on these events, exactly as
    # recompute_published_criterion applies them, and the contrast on THIS outcome
    lab_obs = RPC.classify(real, ao, RPC.strat_nam_150()).astype(float)
    rep_obs = [("Karpechko criterion, CPC AO, days 8-52",
                *matched_contrast(mine, lab_obs))]
    r_obs, pool_obs = analyse("OBSERVATIONS (CPC AO, 43 events)", ao, real, clim,
                              clean_idx, doys, _win_anom, rng, G.draw_clean,
                              reported=rep_obs)
    res["results"]["observations"] = r_obs
    res["validation_known_truth"] = validate_known_truth(
        pool_obs, r_obs["n_events"], rng)

    # ------------------------------------------------ CMIP6
    print("\nloading CMIP6 ensemble ...")
    Yc, Yc_pre, psc, psc_pre, wids, labc = [], [], [], [], [], []
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            m = EP.load_member(f)
            full = m
        except Exception:
            continue
        m = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m) < 2000:
            continue
        on = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
        if len(on) < 15:
            continue
        am = m["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        cln = C6.zone_free_index(am.index, on)
        dy = np.array([t.dayofyear for t in pd.DatetimeIndex(on)])
        Yc.append(C6.anom(am, on, cl, OUT_WIN))
        labc.append(PC.karpechko_surface(am - cl.reindex(am.index.dayofyear).values, on))
        Yc_pre.append(C6.anom(am, on, cl, PRE_WIN))
        wids.append(winter_of(on) + 10000 * len(wids))   # keep winters distinct per member
        for _ in range(20):
            p = C6.draw_pseudo(cln, dy, rng)
            if len(p) >= 10:
                psc.append(C6.anom(am, p, cl, OUT_WIN))
                psc_pre.append(C6.anom(am, p, cl, PRE_WIN))

    if Yc:
        Yall = np.concatenate(Yc)
        Ypre = np.concatenate(Yc_pre)
        pool = np.concatenate(psc)
        pool_pre = np.concatenate(psc_pre)
        wid = np.concatenate(wids)
        print("\n" + "=" * 74)
        print(f"=== CMIP6 ENSEMBLE ===")
        print(f"  events {int(np.isfinite(Yall).sum()):,}, "
              f"pseudo {int(np.isfinite(pool).sum()):,}")
        s2 = sigma_f2(Yall, pool)
        bs, _ = boot_ci(Yall, wid, pool, rng)
        ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        s2p = sigma_f2(Ypre, pool_pre)
        bsp, _ = boot_ci(Ypre, wid, pool_pre, rng)
        cip = [float(np.percentile(bsp, 2.5)), float(np.percentile(bsp, 97.5))]
        vr = float(np.nanvar(Ypre, ddof=1) / np.nanvar(pool_pre, ddof=1))
        print(f"  mean shift S    = {np.nanmean(Yall) - np.nanmean(pool):+.3f} sigma")
        print(f"  Var(SSW) {np.nanvar(Yall, ddof=1):.4f}   "
              f"Var(pseudo) {np.nanvar(pool, ddof=1):.4f}")
        print(f"  s_f^2 = {s2:+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]")
        print(f"  s_f   = {np.sqrt(max(s2, 0)):.3f}  "
              f"(upper 95% {np.sqrt(max(ci[1], 0)):.3f})")
        bsd = boot_paired(Yall, Ypre, wid, pool, pool_pre, rng)
        s2c = s2 - s2p
        cic = [float(np.percentile(bsd, 2.5)), float(np.percentile(bsd, 97.5))]
        print(f"  C1 placebo s_f^2 = {s2p:+.4f} [{cip[0]:+.4f}, {cip[1]:+.4f}] "
              f"-> {'consistent with zero' if cip[0] <= 0 <= cip[1] else 'NON-ZERO confound, corrected below'}")
        print(f"     placebo-CORRECTED s_f^2 = {s2c:+.4f} [{cic[0]:+.4f}, {cic[1]:+.4f}]"
              f"  -> s_f = {np.sqrt(max(s2c, 0)):.3f} "
              f"(upper 95% {np.sqrt(max(cic[1], 0)):.3f})")
        print(f"     {100 * s2p / s2 if s2 else float('nan'):.0f}% of the raw excess "
              f"variance is already present BEFORE onset")
        print(f"  C3 pre-onset variance ratio = {vr:.3f}")
        print(f"\n  === CEILING (uncorrected s_f, at each study's own split q) ===")
        print(f"  {'reported contrast':<50s} {'q':>5s} {'|obs|':>7s} {'ceil95':>7s} "
              f"{'ratio':>7s}  match")
        print("  " + "-" * 84)
        comp = {}
        v, q, n = matched_contrast(Yall, np.concatenate(labc))
        for lab, v, q, n in [("Karpechko conditions 1-2, CMIP6 annular mode", v, q, n)]:
            v = abs(v)
            ch = ceiling(ci[1], q)
            ratio = v / ch if ch > 0 else np.inf
            comp[lab] = {"contrast_abs": round(v, 4), "split_fraction": round(q, 3),
                         "n_events": n,
                         "ceiling_upper95_at_q": round(float(ch), 4),
                         "over_upper_bound_ratio":
                             None if not np.isfinite(ratio) else round(float(ratio), 2),
                         "matched_system": True}
            print(f"  {lab:<50s} {q:5.2f} {v:7.3f} {ch:7.3f} "
                  f"{('inf' if not np.isfinite(ratio) else f'{ratio:.1f}x'):>7s}  YES")
        c_pt, c_hi = ceiling(s2, 0.5), ceiling(ci[1], 0.5)
        print(f"  placebo-CORRECTED ceiling at q=0.5, upper 95% = "
              f"{ceiling(max(cic[1], 0), 0.5):.3f} sigma")
        res["results"]["cmip6"] = {
            "n_events": int(np.isfinite(Yall).sum()),
            "n_pseudo": int(np.isfinite(pool).sum()),
            "mean_shift": round(float(np.nanmean(Yall) - np.nanmean(pool)), 4),
            "var_event": round(float(np.nanvar(Yall, ddof=1)), 4),
            "var_pseudo": round(float(np.nanvar(pool, ddof=1)), 4),
            "sigma_f2": round(s2, 5), "sigma_f2_CI95": [round(c, 5) for c in ci],
            "sigma_f": round(float(np.sqrt(max(s2, 0))), 4),
            "sigma_f_upper95": round(float(np.sqrt(max(ci[1], 0))), 4),
            "placebo_sigma_f2": round(s2p, 5),
            "placebo_CI95": [round(c, 5) for c in cip],
            "placebo_clean": bool(cip[0] <= 0 <= cip[1]),
            "placebo_share_of_raw_pct": round(float(100 * s2p / s2), 1) if s2 else None,
            "corrected_sigma_f2": round(float(s2c), 5),
            "corrected_CI95": [round(c, 5) for c in cic],
            "corrected_sigma_f": round(float(np.sqrt(max(s2c, 0))), 4),
            "corrected_sigma_f_upper95": round(float(np.sqrt(max(cic[1], 0))), 4),
            "corrected_ceiling_upper95": round(float(ceiling(max(cic[1], 0))), 4),
            "pre_onset_variance_ratio": round(vr, 4),
            "ceiling_point": round(float(c_pt), 4),
            "ceiling_upper95": round(float(c_hi), 4),
            "reported_vs_ceiling": comp}

    (RESULTS / "forced_variance_ceiling.json").write_text(
        json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> forced_variance_ceiling.json")


if __name__ == "__main__":
    main()
