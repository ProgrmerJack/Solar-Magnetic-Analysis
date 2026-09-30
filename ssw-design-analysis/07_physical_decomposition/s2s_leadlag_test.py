#!/usr/bin/env python3
"""
s2s_leadlag_test.py
===================
WHERE IN THE FORECAST CHAIN IS THE SHIFT LOST: THE LOWER STRATOSPHERE, OR BELOW IT?

Plan approved 2026-09-30. This design was committed BEFORE any forecast of
100 hPa geopotential height or of 10 hPa wind beyond lead 15 was retrieved.

QUESTION
  Operational ensembles started 2-9 days before SSWs place the observed days
  +8..+25 polar-cap surface outcome too far on the downward side (P', mean rank
  0.40 against 0.54). Is that because the forecast lower-stratospheric anomaly is
  too weak or too short-lived (a stratospheric persistence deficit, as reported
  for most S2S systems by Garfinkel et al. 2025), or because, given the forecast
  lower stratosphere, the surface responds too little (a coupling deficit)?

DATA
  Forecasts: the ten P' systems, model versions, starts and members; polar-cap
  (60-90N, cos-lat, the 60N row included) geopotential height at 100 hPa, 00 UTC,
  leads 1-34 (acquire_s2s_reforecasts.py --var gh100); zonal-mean u at 10 hPa, 60N,
  leads 16-34 (--var u10_long), joined to the cached leads 1-15. Surface outcome as
  in P'. Observations: ERA5 (WeatherBench 2) 100 hPa geopotential height, 00 UTC,
  the same cap (acquire_era5_psl_cap.py --z100); ERA5 polar-cap msl as in P'.
  Anomalies: leave-one-year-out, as in P'. Z100 anomalies are sign-reversed
  (B = -anomaly), so that NEGATIVE is the weak-vortex direction, like A.

TESTS (starts 2-9 d before onset; the 17 P' events; multi-model mean of the nine
confirmatory systems; the P' calendar-window null with the quorum rule;
one-sided p = P(null <= observed))
  L1 (primary)  rank of the observed days +8..+25 mean B within the ensemble: is
                the lower-stratospheric anomaly after SSWs under-forecast?
  L2 (primary)  conditional surface rank: within each ensemble, regress members' A
                (surface, days +8..+25) on their B (same window); the conditional
                forecast at the OBSERVED B is a + b*B_obs plus the members'
                residuals; rank of the observed A in it. Is the surface
                under-forecast even given the lower stratosphere?
  L3 (descriptive) decomposition of the ensemble-mean surface error per event:
                A_obs - A_ens = b (B_obs - B_ens) + remainder; multi-model means of
                both parts with 10,000-resample event-bootstrap intervals.
  L4 (secondary) the same as L1 for u(10 hPa, 60N) over days +8..+25.
  Holm correction over L1 and L2 (this diagnostic family only; they are not added
  to the paper's primary family, whose register is unchanged).
  Reading, fixed now:
    L1 low, L2 ~ null   -> the shift is lost in the stratosphere (persistence);
    L1 ~ null, L2 low   -> it is lost below 100 hPa (coupling);
    both low            -> both;  neither -> the diagnostic does not locate it.
  ECMWF (discovery system) is reported separately and is not in the primary mean.

IMPLEMENTATION NOTES (fixed before any 100 hPa forecast was read)
  Surface A and B are both days +8..+25 means (leads k+8..k+25). The within-ensemble
  regression uses all members of the start (4-11); a start whose members' B has
  zero variance is skipped for L2/L3. An event's value per system is the mean over
  its starts; multi-model mean and null exactly as P' (quorum of half). L4's
  observation is the NCEP-NCAR daily-mean zonal-mean u at 10 hPa, 60N
  (extend_ncep_presatellite cache, as validate_ssw_detector.py), against 00 UTC
  forecasts; forecasts join the cached leads 1-15 with leads 16-34.

Output: results/current/6_predictability/s2s_leadlag_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ING))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_leadlag_test"
N_NULL = 10000
N_BOOT = 10000
WIN = (8, 25)
OBS_Z = ING / "era5_z100_cap_00utc.parquet"


def tag(c):
    return "ecmf" if c == "ecmwf" else c


def anom(fc, col):
    """Leave-one-year-out forecast anomaly at (start key, lead), as T.load."""
    fc = fc.copy()
    fc["init"] = pd.to_datetime(fc["init"]); fc["hyear"] = fc["init"].dt.year
    fc["md"] = pd.to_datetime(fc["model_date"]).dt.strftime("%m-%d")
    g = fc.groupby(["md", "lead_day", "hyear"])[col].agg(["sum", "count"])
    tot = g.groupby(level=[0, 1]).sum()
    loo = tot.reindex(g.index.droplevel(2)).values - g.values
    clim = pd.Series(loo[:, 0] / loo[:, 1], index=g.index)
    fc["anom"] = fc[col].values - clim.reindex(pd.MultiIndex.from_frame(fc[["md", "lead_day", "hyear"]])).values
    return fc


def obs_series():
    ob = pd.read_parquet(T.OBS); ob["time"] = pd.to_datetime(ob["time"])
    psl = ob[ob["time"].dt.hour == 0].set_index("time")["psl_cap_N"]
    z = pd.read_parquet(OBS_Z); z = z.set_index(pd.to_datetime(z["time"]))["z100_cap_N"]
    import extend_ncep_presatellite as EN
    u = pd.concat([EN.series("uwnd", y, 10) for y in range(1998, 2024)])["uwnd_ms_60N_10hPa"]
    u.index = pd.to_datetime(u.index).tz_localize(None).normalize()
    return psl, z, u[~u.index.duplicated()].sort_index()


class Sys:
    def __init__(self, c, psl_ob, z_ob, u_ob):
        a = anom(pd.read_parquet(MM.FILES[c]), "psl_cap_N")
        b = anom(pd.read_parquet(ING / f"s2s_{tag(c)}_z100_cap.parquet"), "z100_cap_N")
        u = pd.concat([pd.read_parquet(ING / f"s2s_{tag(c)}_u10_60N.parquet"),
                       pd.read_parquet(ING / f"s2s_{tag(c)}_u10_60N_long.parquet")], ignore_index=True)
        u = anom(u, "u10_60N")
        self.hyears = sorted(int(y) for y in a["hyear"].unique())
        self.min_other = max(8, int(0.75 * (len(self.hyears) - 1)))
        self.A = {pd.Timestamp(i): g.pivot(index="member", columns="lead_day", values="anom") for i, g in a.groupby("init")}
        self.B = {pd.Timestamp(i): g.pivot(index="member", columns="lead_day", values="anom") for i, g in b.groupby("init")}
        self.U = {pd.Timestamp(i): g.pivot(index="member", columns="lead_day", values="anom") for i, g in u.groupby("init")}
        self.inits = sorted(set(self.A) & set(self.B))
        self.ob = {"A": psl_ob, "B": z_ob, "U": u_ob}
        self._c = {}

    def _obs(self, key, onset):
        o = T.obs_anom(self.ob[key], [onset], range(WIN[0], WIN[1] + 1), self.hyears, self.min_other)[onset]
        return None if np.isnan(o).any() else float(o.mean())

    def at(self, p):
        p = pd.Timestamp(p)
        if p in self._c:
            return self._c[p]
        rows = []
        for i in self.inits:
            k = (p - i).days
            if not (MM.K_SHORT[0] <= k <= MM.K_SHORT[1]):
                continue
            leads = list(range(k + WIN[0], k + WIN[1] + 1))
            fa, fb = self.A[i], self.B[i]
            if not set(leads) <= set(fa.columns) or not set(leads) <= set(fb.columns):
                continue
            Am, Bm = -fa[leads].mean(axis=1), -fb[leads].mean(axis=1)
            j = Am.index.intersection(Bm.index)
            Am, Bm = Am[j].values, Bm[j].values
            if len(Am) < 3 or np.isnan(Am).any() or np.isnan(Bm).any():
                continue
            Ao, Bo = self._obs("A", p), self._obs("B", p)
            if Ao is None or Bo is None:
                continue
            Ao, Bo = -Ao, -Bo
            n = len(Am)
            pit = lambda m, o: float((np.sum(m < o) + 0.5 * np.sum(m == o) + 0.5) / (n + 1))
            r = {"pit_B": pit(Bm, Bo), "A_err": Ao - Am.mean(), "B_err": Bo - Bm.mean()}
            if Bm.std() > 0:
                bcoef = np.polyfit(Bm, Am, 1)
                res = Am - np.polyval(bcoef, Bm)
                r.update(pit_cond=pit(np.polyval(bcoef, Bo) + res, Ao), b=float(bcoef[0]),
                         strat_part=float(bcoef[0] * (Bo - Bm.mean())))
                r["remainder"] = r["A_err"] - r["strat_part"]
            fu = self.U.get(i)
            if fu is not None and set(leads) <= set(fu.columns):
                Um = -fu[leads].mean(axis=1).values
                Uo = self._obs("U", p)
                if Uo is not None and not np.isnan(Um).any():
                    r["pit_U"] = pit(Um, -Uo)
            rows.append(r)
        out = None
        if rows:
            out = {k: float(np.mean([r[k] for r in rows if k in r])) for k in
                   ("pit_B", "pit_cond", "A_err", "B_err", "strat_part", "remainder", "pit_U", "b")
                   if any(k in r for r in rows)}
            out["n_starts"] = len(rows)
        self._c[p] = out
        return out


def main():
    rng = np.random.default_rng(zlib.crc32(NAME.encode()))
    psl, z, u = obs_series()
    cat = load_catalogue("primary")
    cen = {}
    for c in [MM.DISCOVERY] + MM.CONFIRM:
        try:
            cen[c] = Sys(c, psl, z, u)
        except FileNotFoundError as e:
            print(f"  {c}: missing {e.filename}")
    res = {"plan_approved": "2026-09-30", "registered_commit": "4c8b3f2", "n_null": N_NULL,
           "n_boot": N_BOOT, "seed": zlib.crc32(NAME.encode()), "window": list(WIN),
           "systems": list(cen), "tests": {}}

    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))

    for label, group in (("confirmatory", [c for c in MM.CONFIRM if c in cen]),
                         ("ecmwf_discovery", [c for c in [MM.DISCOVERY] if c in cen])):
        if not group:
            continue
        cover = {o: [c for c in group if cen[c].at(o) is not None] for o in cat}
        ev = [o for o in cat if cover[o]]
        years = sorted({y for c in group for y in cen[c].hyears} | {y + 1 for c in group for y in cen[c].hyears})
        need = {o: int(np.ceil(len(cover[o]) / 2)) for o in ev}
        pools = {o: [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
                     for p in [o + pd.DateOffset(years=y - o.year) + pd.Timedelta(days=dd)]
                     if zone_free(p) and sum(cen[c].at(p) is not None for c in cover[o]) >= need[o]] for o in ev}

        def mm(p, systems, key):
            v = [cen[c].at(p)[key] for c in systems if cen[c].at(p) is not None and key in cen[c].at(p)]
            return float(np.mean(v)) if v else np.nan

        out = {"n_events": len(ev)}
        for key, nm in (("pit_B", "L1_z100_rank"), ("pit_cond", "L2_conditional_surface_rank"),
                        ("pit_U", "L4_u10_rank")):
            obs = np.nanmean([mm(o, cover[o], key) for o in ev])
            null = np.array([np.nanmean([mm(pools[o][rng.integers(len(pools[o]))], cover[o], key) for o in ev])
                             for _ in range(N_NULL)])
            out[nm] = {"mean_rank": round(float(obs), 4), "null_mean": round(float(np.nanmean(null)), 4),
                       "null_q025_q975": [round(float(q), 4) for q in np.nanquantile(null, [0.025, 0.975])],
                       "p": round(float(np.nanmean(null <= obs)), 4)}
        # L3 decomposition, event bootstrap
        E = np.array([[mm(o, cover[o], k) for k in ("A_err", "strat_part", "remainder")] for o in ev])
        E = E[np.isfinite(E).all(1)]
        bs = np.array([E[rng.integers(0, len(E), len(E))].mean(0) for _ in range(N_BOOT)])
        out["L3_decomposition_Pa"] = {k: {"mean": round(float(E[:, j].mean()), 2),
                                          "ci95": [round(float(q), 2) for q in np.percentile(bs[:, j], [2.5, 97.5])]}
                                      for j, k in enumerate(("surface_error", "via_100hPa", "remainder"))}
        out["L3_n_events"] = int(len(E))
        out["mean_within_ensemble_slope_Pa_per_gpm"] = round(float(np.nanmean([mm(o, cover[o], "b") for o in ev])), 4)
        if label == "confirmatory":
            raw = {"L1": out["L1_z100_rank"]["p"], "L2": out["L2_conditional_surface_rank"]["p"]}
            (k1, p1), (k2, p2) = sorted(raw.items(), key=lambda kv: kv[1])
            a1 = min(1.0, 2 * p1)
            out["holm_L1_L2"] = {k1: round(a1, 4), k2: round(max(a1, p2), 4)}
        res["tests"][label] = out
        print(label, json.dumps(out), flush=True)
    (RESULTS / "s2s_leadlag_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_leadlag_test.json")


if __name__ == "__main__":
    main()
