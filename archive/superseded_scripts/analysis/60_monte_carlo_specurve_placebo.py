"""
60_monte_carlo_specurve_placebo.py
Nature Geoscience manuscript: SSW → avalanche hazard

Three analyses:
  1. Monte Carlo specification-curve significance (10k iterations)
  2. Placebo / falsification windows
  3. Leave-two-out event influence

Outputs → data/results/r55_specurve_placebo.json
"""

import json, warnings, itertools, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
np.random.seed(42)

ROOT = Path(__file__).resolve().parents[2]
PANEL_PATH = ROOT / "data" / "processed" / "analysis_panel_v2.parquet"
SSW_PATH = ROOT / "data" / "processed" / "atmospheric" / "ssw_catalog.parquet"
OUT_PATH = ROOT / "data" / "results" / "r55_specurve_placebo.json"

N_MC = 10_000
WINDOWS = [5, 10, 15, 20, 25, 30]
DOY_BWS = [1, 3, 5, 7, 10]
AVAL_COL = "dry_natural_size_1234"

# ── Load & pre-compute numpy arrays for speed ───────────────────────────

def load_data():
    panel = pd.read_parquet(PANEL_PATH)
    if panel.index.tz is not None:
        panel.index = panel.index.tz_localize(None)
    panel = panel[panel["is_winter"] == 1].copy()

    ssw_cat = pd.read_parquet(SSW_PATH)
    ssw_dates = ssw_cat.index
    if ssw_dates.tz is not None:
        ssw_dates = ssw_dates.tz_localize(None)
    ssw_dates = ssw_dates[(ssw_dates >= panel.index.min()) & (ssw_dates <= panel.index.max())]
    return panel, ssw_dates


print("Loading data …")
t0 = time.time()
panel, ssw_dates = load_data()
n_events = len(ssw_dates)
all_winters = sorted(panel["winter_id"].unique())
n_winters = len(all_winters)
winter_to_idx = {w: i for i, w in enumerate(all_winters)}

# Convert to fast numpy arrays
dates_np = panel.index.values.astype("datetime64[D]").astype(np.int64)  # days since epoch
aval_np = panel[AVAL_COL].fillna(0).values.astype(np.float64)
doy_np = panel["day_of_year"].values.astype(np.int32)
wid_np = panel["winter_id"].map(winter_to_idx).values.astype(np.int32)
u10_np = panel["ncep_u_10hpa"].values.astype(np.float64)

# Date-to-index lookup
date_to_pos = {d: i for i, d in enumerate(dates_np)}

# Indices per winter
winter_indices = {i: np.where(wid_np == i)[0] for i in range(n_winters)}

ssw_dates_list = list(ssw_dates)
ssw_days_epoch = np.array([np.datetime64(d, "D").astype(np.int64) for d in ssw_dates_list])
ssw_doys = np.array([pd.Timestamp(d).day_of_year for d in ssw_dates_list])
ssw_winter_idx = np.array([
    wid_np[np.argmin(np.abs(dates_np - sd))] for sd in ssw_days_epoch
])

print(f"  Panel: {len(panel)} winter days, {n_winters} winters")
print(f"  SSW events: {n_events}")
print(f"  Pre-processing took {time.time()-t0:.1f}s")


# ── Vectorized RR computation ───────────────────────────────────────────

def compute_rr_batch(onset_epochs, onset_doys, onset_wids, window, doy_bw):
    """Compute per-event RR for a batch of onset dates. All numpy."""
    n = len(onset_epochs)
    rr = np.full(n, np.nan)
    for i in range(n):
        d_epoch = onset_epochs[i]
        doy = onset_doys[i]
        w_idx = onset_wids[i]

        # Observed: days within ±window of onset
        lo = d_epoch - window
        hi = d_epoch + window
        obs_mask = (dates_np >= lo) & (dates_np <= hi)
        observed = aval_np[obs_mask].sum()

        # Expected: DOY-matched days in non-event winters
        doy_lo = doy - doy_bw
        doy_hi = doy + doy_bw
        if doy_lo < 1:
            doy_match = (doy_np >= (doy_lo + 366)) | (doy_np <= doy_hi)
        elif doy_hi > 366:
            doy_match = (doy_np >= doy_lo) | (doy_np <= (doy_hi - 366))
        else:
            doy_match = (doy_np >= doy_lo) & (doy_np <= doy_hi)

        ctrl_mask = doy_match & (wid_np != w_idx)
        ctrl_vals = aval_np[ctrl_mask]
        n_ctrl = len(ctrl_vals)
        if n_ctrl == 0:
            continue

        ctrl_daily_rate = ctrl_vals.sum() / n_ctrl
        expected = ctrl_daily_rate * (2 * window + 1)
        if expected > 0:
            rr[i] = observed / expected
    return rr


def geometric_mean_rr(rr_arr):
    valid = rr_arr[np.isfinite(rr_arr) & (rr_arr > 0)]
    if len(valid) == 0:
        return np.nan
    return float(np.exp(np.mean(np.log(valid))))


def sign_count(rr_arr):
    valid = rr_arr[np.isfinite(rr_arr)]
    return int(np.sum(valid < 1))


# ── Actual SSW results ──────────────────────────────────────────────────

actual_rr = compute_rr_batch(ssw_days_epoch, ssw_doys, ssw_winter_idx, 15, 3)
actual_gmrr = geometric_mean_rr(actual_rr)
actual_n_dec = sign_count(actual_rr)
actual_sign_p = float(stats.binomtest(actual_n_dec, int(np.sum(np.isfinite(actual_rr))), 0.5).pvalue)
print(f"  Actual: gmRR={actual_gmrr:.4f}, n_decrease={actual_n_dec}/{n_events}, P={actual_sign_p:.4f}")


# ═════════════════════════════════════════════════════════════════════════
#  ANALYSIS 1 — Monte Carlo specification-curve significance
# ═════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"ANALYSIS 1: Monte Carlo spec-curve ({N_MC} iterations)")
print("=" * 60)
t1 = time.time()

mc_n_decrease = np.zeros(N_MC, dtype=int)
mc_gmrr = np.zeros(N_MC, dtype=float)
mc_specs_lt1_count = np.zeros(N_MC, dtype=int)
n_specs = len(WINDOWS) * len(DOY_BWS)

for i in range(N_MC):
    if i % 2000 == 0:
        elapsed = time.time() - t1
        eta = (elapsed / max(i, 1)) * (N_MC - i)
        print(f"  MC {i}/{N_MC}  ({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    # Draw 16 random winter dates
    rand_w = np.random.randint(0, n_winters, size=n_events)
    rand_pos = np.array([np.random.choice(winter_indices[w]) for w in rand_w])
    rand_epochs = dates_np[rand_pos]
    rand_doys = doy_np[rand_pos]

    # Primary spec: window=15, doy_bw=3
    rr = compute_rr_batch(rand_epochs, rand_doys, rand_w, 15, 3)
    mc_n_decrease[i] = sign_count(rr)
    mc_gmrr[i] = geometric_mean_rr(rr)

    # All 30 specs
    s_lt1 = 0
    for w in WINDOWS:
        for bw in DOY_BWS:
            rr_s = compute_rr_batch(rand_epochs, rand_doys, rand_w, w, bw)
            gm = geometric_mean_rr(rr_s)
            if np.isfinite(gm) and gm < 1:
                s_lt1 += 1
    mc_specs_lt1_count[i] = s_lt1

mc_frac = mc_specs_lt1_count / n_specs
mc_all_lt1 = mc_specs_lt1_count == n_specs

p_14_of_16 = float(np.mean(mc_n_decrease >= 14))
p_16_of_16 = float(np.mean(mc_n_decrease >= 16))
p_actual_dec = float(np.mean(mc_n_decrease >= actual_n_dec))
p_100_spec = float(np.mean(mc_all_lt1))
p_90_spec = float(np.mean(mc_frac >= 0.90))

print(f"\n  Completed in {time.time()-t1:.0f}s")
print(f"  P(≥{actual_n_dec}/{n_events} decrease) = {p_actual_dec:.6f}")
print(f"  P(≥14/16 decrease)       = {p_14_of_16:.6f}")
print(f"  P(16/16 decrease)        = {p_16_of_16:.6f}")
print(f"  P(100% spec consistency) = {p_100_spec:.6f}")
print(f"  P(≥90% spec consistency) = {p_90_spec:.6f}")
print(f"  Null gmRR: mean={np.nanmean(mc_gmrr):.4f}, sd={np.nanstd(mc_gmrr):.4f}")

mc_result = {
    "n_iterations": int(N_MC),
    "n_specs_per_iteration": int(n_specs),
    "actual_gmrr": round(float(actual_gmrr), 6),
    "actual_n_decrease": int(actual_n_dec),
    "actual_sign_p": round(actual_sign_p, 6),
    "p_gte_actual_decrease": round(p_actual_dec, 6),
    "p_14_of_16_sign": round(p_14_of_16, 6),
    "p_16_of_16_sign": round(p_16_of_16, 6),
    "p_100pct_spec_consistency": round(p_100_spec, 6),
    "p_90pct_spec_consistency": round(p_90_spec, 6),
    "null_distribution_summary": {
        "gmrr_mean": round(float(np.nanmean(mc_gmrr)), 6),
        "gmrr_median": round(float(np.nanmedian(mc_gmrr)), 6),
        "gmrr_sd": round(float(np.nanstd(mc_gmrr)), 6),
        "gmrr_5th": round(float(np.nanpercentile(mc_gmrr, 5)), 6),
        "gmrr_95th": round(float(np.nanpercentile(mc_gmrr, 95)), 6),
        "n_decrease_mean": round(float(mc_n_decrease.mean()), 4),
        "n_decrease_sd": round(float(mc_n_decrease.std()), 4),
        "spec_frac_lt1_mean": round(float(mc_frac.mean()), 6),
        "spec_frac_lt1_sd": round(float(mc_frac.std()), 6),
    },
}


# ═════════════════════════════════════════════════════════════════════════
#  ANALYSIS 2 — Placebo / Falsification
# ═════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("ANALYSIS 2: Placebo / falsification windows")
print("=" * 60)

def placebo_summary(rr_arr, label):
    valid = rr_arr[np.isfinite(rr_arr)]
    n_test = len(valid)
    nd = int(np.sum(valid < 1))
    gm = geometric_mean_rr(rr_arr)
    sp = float(stats.binomtest(nd, n_test, 0.5).pvalue) if n_test > 0 else 1.0
    print(f"    {label}: gmRR={gm:.4f}, n_decrease={nd}/{n_test}, P={sp:.4f}")
    return {
        "mean_RR": round(float(gm), 6) if np.isfinite(gm) else None,
        "n_tested": n_test,
        "n_decrease": nd,
        "sign_P": round(sp, 6),
    }


def epochs_doys_wids(date_list):
    epochs = np.array([np.datetime64(pd.Timestamp(d), "D").astype(np.int64) for d in date_list])
    doys = np.array([pd.Timestamp(d).day_of_year for d in date_list])
    wids = np.array([
        wid_np[np.argmin(np.abs(dates_np - e))] if np.min(np.abs(dates_np - e)) < 60 else 0
        for e in epochs
    ])
    return epochs, doys, wids


# 2a: +60 day shift
print("  2a: +60 day shift")
shifted_p60 = [d + pd.Timedelta(days=60) for d in ssw_dates_list]
ep, dy, wi = epochs_doys_wids(shifted_p60)
rr_p60 = compute_rr_batch(ep, dy, wi, 15, 3)
res_p60 = placebo_summary(rr_p60, "+60d")

# 2b: -60 day shift
print("  2b: -60 day shift")
shifted_m60 = [d - pd.Timedelta(days=60) for d in ssw_dates_list]
ep, dy, wi = epochs_doys_wids(shifted_m60)
rr_m60 = compute_rr_batch(ep, dy, wi, 15, 3)
res_m60 = placebo_summary(rr_m60, "-60d")

# 2c: Random dates (1000 iterations)
print("  2c: Random date placebo (1000 iters)")
N_RAND = 1000
rand_frac_dec = np.zeros(N_RAND)
for j in range(N_RAND):
    rw = np.random.randint(0, n_winters, size=n_events)
    rp = np.array([np.random.choice(winter_indices[w]) for w in rw])
    rr_r = compute_rr_batch(dates_np[rp], doy_np[rp], rw, 15, 3)
    vr = rr_r[np.isfinite(rr_r)]
    rand_frac_dec[j] = np.sum(vr < 1) / max(len(vr), 1)

rand_p_gte14 = float(np.mean(rand_frac_dec >= 14 / 16))
print(f"    Mean frac decrease = {rand_frac_dec.mean():.4f}, P(≥14/16) = {rand_p_gte14:.6f}")

res_rand = {
    "n_iterations": int(N_RAND),
    "mean_fraction_decrease": round(float(rand_frac_dec.mean()), 6),
    "sd_fraction_decrease": round(float(rand_frac_dec.std()), 6),
    "p_gte_14_of_16": round(rand_p_gte14, 6),
}

# 2d: Non-SSW winters (same DOY as SSW, but in winters without SSW)
print("  2d: Non-SSW winters")
ssw_winter_set = set(ssw_winter_idx.tolist())
non_ssw_winter_list = [i for i in range(n_winters) if i not in ssw_winter_set]
print(f"    SSW winters: {len(ssw_winter_set)}, non-SSW: {len(non_ssw_winter_list)}")

rr_non_ssw_all = []
for si in range(n_events):
    doy_target = ssw_doys[si]
    for nw in non_ssw_winter_list:
        idxs = winter_indices[nw]
        doy_match = np.where(doy_np[idxs] == doy_target)[0]
        if len(doy_match) > 0:
            pos = idxs[doy_match[0]]
            rr_val = compute_rr_batch(
                dates_np[pos:pos+1], doy_np[pos:pos+1],
                np.array([nw]), 15, 3
            )[0]
            if np.isfinite(rr_val):
                rr_non_ssw_all.append(rr_val)

rr_non_ssw_arr = np.array(rr_non_ssw_all) if rr_non_ssw_all else np.array([np.nan])
res_non = placebo_summary(rr_non_ssw_arr, "non-SSW winters")

# 2e: Weak vortex non-SSW
print("  2e: Weak vortex (non-SSW)")
p20 = float(np.nanpercentile(u10_np[np.isfinite(u10_np)], 20))
weak_mask = np.isfinite(u10_np) & (u10_np > 0) & (u10_np < p20)

# Exclude ±30d around actual SSW events
for sd in ssw_days_epoch:
    vicinity = (dates_np >= sd - 30) & (dates_np <= sd + 30)
    weak_mask = weak_mask & ~vicinity

weak_positions = np.where(weak_mask)[0]
print(f"    Weak-vortex days: {len(weak_positions)} (u10 threshold={p20:.2f})")

if len(weak_positions) >= n_events:
    np.random.shuffle(weak_positions)
    sel = weak_positions[:n_events]
    rr_weak = compute_rr_batch(dates_np[sel], doy_np[sel], wid_np[sel], 15, 3)
else:
    rr_weak = np.array([np.nan])

res_weak = placebo_summary(rr_weak, "weak vortex")
res_weak["u10_threshold"] = round(p20, 4)

placebo_result = {
    "actual_ssw": {
        "mean_RR": round(float(actual_gmrr), 6),
        "n_decrease": int(actual_n_dec),
        "sign_P": round(actual_sign_p, 6),
    },
    "shifted_plus60": res_p60,
    "shifted_minus60": res_m60,
    "random_dates": res_rand,
    "non_ssw_winters": res_non,
    "weak_vortex_non_ssw": res_weak,
}


# ═════════════════════════════════════════════════════════════════════════
#  ANALYSIS 3 — Leave-Two-Out Event Influence
# ═════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("ANALYSIS 3: Leave-two-out event influence")
print("=" * 60)

combos = list(itertools.combinations(range(n_events), 2))
n_combos = len(combos)
print(f"  C({n_events},2) = {n_combos} combinations")

l2o_gmrr = np.zeros(n_combos)
l2o_n_dec = np.zeros(n_combos, dtype=int)
l2o_sign_p = np.zeros(n_combos)

for ci, (a, b) in enumerate(combos):
    keep = [k for k in range(n_events) if k != a and k != b]
    rr = compute_rr_batch(ssw_days_epoch[keep], ssw_doys[keep], ssw_winter_idx[keep], 15, 3)
    valid = rr[np.isfinite(rr)]
    l2o_gmrr[ci] = geometric_mean_rr(rr)
    nd = int(np.sum(valid < 1))
    l2o_n_dec[ci] = nd
    l2o_sign_p[ci] = float(stats.binomtest(nd, len(valid), 0.5).pvalue) if len(valid) > 0 else 1.0

worst_idx = int(np.argmax(l2o_sign_p))
worst_combo = combos[worst_idx]
min_dec_idx = int(np.argmin(l2o_n_dec))
min_dec_combo = combos[min_dec_idx]

all_gmrr_lt1 = bool(np.all(l2o_gmrr[np.isfinite(l2o_gmrr)] < 1))
all_sig_010 = bool(np.all(l2o_sign_p < 0.10))
all_sig_005 = bool(np.all(l2o_sign_p < 0.05))

print(f"  Min sign P   = {l2o_sign_p.min():.6f}")
print(f"  Max sign P   = {l2o_sign_p.max():.6f}")
print(f"  Min n_dec    = {l2o_n_dec.min()}")
print(f"  All gmRR<1?  = {all_gmrr_lt1}")
print(f"  All P<0.05?  = {all_sig_005}")
print(f"  Worst drop   = {str(ssw_dates_list[worst_combo[0]])[:10]}, {str(ssw_dates_list[worst_combo[1]])[:10]}")

l2o_result = {
    "n_combinations": int(n_combos),
    "n_events_total": int(n_events),
    "min_sign_p": round(float(l2o_sign_p.min()), 6),
    "max_sign_p": round(float(l2o_sign_p.max()), 6),
    "median_sign_p": round(float(np.median(l2o_sign_p)), 6),
    "min_n_decrease": int(l2o_n_dec.min()),
    "max_n_decrease": int(l2o_n_dec.max()),
    "min_gmrr": round(float(np.nanmin(l2o_gmrr)), 6),
    "max_gmrr": round(float(np.nanmax(l2o_gmrr)), 6),
    "all_gmrr_lt1": all_gmrr_lt1,
    "all_remain_significant_p010": all_sig_010,
    "all_remain_significant_p005": all_sig_005,
    "worst_case": {
        "dropped_events": [
            str(ssw_dates_list[worst_combo[0]])[:10],
            str(ssw_dates_list[worst_combo[1]])[:10],
        ],
        "sign_p": round(float(l2o_sign_p[worst_idx]), 6),
        "gmrr": round(float(l2o_gmrr[worst_idx]), 6),
        "n_decrease": int(l2o_n_dec[worst_idx]),
        "n_tested": int(n_events - 2),
    },
    "min_decrease_case": {
        "dropped_events": [
            str(ssw_dates_list[min_dec_combo[0]])[:10],
            str(ssw_dates_list[min_dec_combo[1]])[:10],
        ],
        "n_decrease": int(l2o_n_dec[min_dec_idx]),
        "gmrr": round(float(l2o_gmrr[min_dec_idx]), 6),
        "sign_p": round(float(l2o_sign_p[min_dec_idx]), 6),
    },
}


# ═════════════════════════════════════════════════════════════════════════
#  SAVE
# ═════════════════════════════════════════════════════════════════════════
results = {
    "monte_carlo_spec_curve": mc_result,
    "placebo_falsification": placebo_result,
    "leave_two_out": l2o_result,
}

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*60}")
print(f"All results saved → {OUT_PATH.name}")
print(f"Total time: {time.time()-t0:.0f}s")
print("=" * 60)

