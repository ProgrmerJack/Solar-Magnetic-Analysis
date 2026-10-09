#!/usr/bin/env python3
"""
acquire_ghcnd_neurasia.py
=========================
Input of 07_physical_decomposition/early_transfer_retest.py (registered 2026-10-09, committed
before any station value was read): a station-based northern-Eurasian (50-65 N, 10-130 E)
daily-mean temperature anomaly from GHCN-Daily (Menne et al. 2012), NOAA NCEI, anonymous:
  https://www.ncei.noaa.gov/pub/data/ghcn/daily/{ghcnd-inventory.txt, ghcnd-stations.txt,
  by_station/<ID>.csv.gz}
Stations: inside the box, with TMAX and TMIN whose inventory first year is <= 1941 and last
year >= 2022 (fixed before any value was read). Daily mean = (TMAX + TMIN) / 2, both present
with a blank quality flag. Per station and per base period, the anomaly from that station's
31-day-smoothed day-of-year mean (base 1941-1958 for the early tests; base 1959-2022 for the
validation gate); a station-base needs >= 70% of its base-period November-March days.
Regional series: station anomalies averaged within 5 deg x 10 deg cells, then cos(latitude of
the cell centre)-weighted over cells; a day is kept if >= 70% of the cells with a station
report. Writes ghcnd_neurasia_daily.parquet (date, anom_b1941_1958, anom_b1959_2022,
n_stations, n_cells) and the station list ghcnd_neurasia_stations.csv; raw station files in
raw/ghcnd/.
"""
import gzip
import io
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "ghcnd"
URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/"
OUT = HERE / "ghcnd_neurasia_daily.parquet"
LST = HERE / "ghcnd_neurasia_stations.csv"
BOX = (50, 65, 10, 130)
BASES = {"anom_b1941_1958": (1941, 1958), "anom_b1959_2022": (1959, 2022)}
NDJFM = (11, 12, 1, 2, 3)


def fetch(name):
    p = RAW / name
    if not p.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL + ("by_station/" if name.endswith(".csv.gz") else "") + name, timeout=300) as r:
            p.write_bytes(r.read())
    return p


def stations():
    inv = pd.read_fwf(fetch("ghcnd-inventory.txt"), colspecs=[(0, 11), (12, 20), (21, 30), (31, 35), (36, 40), (41, 45)],
                      names=["id", "lat", "lon", "el", "first", "last"])
    x = inv[inv.lat.between(BOX[0], BOX[1]) & inv.lon.between(BOX[2], BOX[3])]
    tx = x[x.el == "TMAX"].set_index("id"); tn = x[x.el == "TMIN"].set_index("id")
    j = tx.join(tn[["first", "last"]], rsuffix="_n", how="inner")
    s = j[(j["first"] <= 1941) & (j["first_n"] <= 1941) & (j["last"] >= 2022) & (j["last_n"] >= 2022)]
    return s[["lat", "lon"]].reset_index()


def station_series(sid):
    raw = gzip.decompress(fetch(f"{sid}.csv.gz").read_bytes()).decode()
    d = pd.read_csv(io.StringIO(raw), header=None, usecols=[1, 2, 3, 5], names=["date", "el", "v", "q"], dtype={"q": str},
                    keep_default_na=False)
    d = d[d.el.isin(["TMAX", "TMIN"]) & (d.q == "")]
    p = d.pivot_table(index="date", columns="el", values="v", aggfunc="first")
    if not {"TMAX", "TMIN"} <= set(p.columns):
        return None
    t = ((p.TMAX + p.TMIN) / 20.0).dropna()
    t.index = pd.to_datetime(t.index.astype(str), format="%Y%m%d")
    return t


def anomaly(t, base):
    b = t[(t.index.year >= base[0]) & (t.index.year <= base[1])]
    need = sum(1 for d in pd.date_range(f"{base[0]}-01-01", f"{base[1]}-12-31") if d.month in NDJFM)
    if (b.index.month.isin(NDJFM)).sum() < 0.7 * need:
        return None
    doy = b.groupby(b.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    return t - doy.reindex(t.index.dayofyear).values


def main():
    st = stations()
    st["cell"] = list(zip((st.lat // 5).astype(int), (st.lon // 10).astype(int)))
    print(f"{len(st)} stations, {st.cell.nunique()} cells", flush=True)
    with ThreadPoolExecutor(8) as ex:
        ser = dict(zip(st.id, ex.map(station_series, st.id)))
    out = {}
    for col, base in BASES.items():
        an = {s: anomaly(t, base) for s, t in ser.items() if t is not None}
        an = {s: a for s, a in an.items() if a is not None}
        df = pd.DataFrame(an)
        cells = st.set_index("id").loc[df.columns, "cell"]
        cm = df.T.groupby(cells.values).mean().T                  # cell means (NaN-skipping)
        w = np.cos(np.radians([5 * c[0] + 2.5 for c in cm.columns]))
        ok = cm.notna()
        reg = (cm.fillna(0) * w).sum(axis=1) / (ok * w).sum(axis=1)
        reg[ok.sum(axis=1) < 0.7 * cm.shape[1]] = np.nan
        out[col] = reg
        out[f"n_stations_{col[5:]}"] = df.notna().sum(axis=1)
        out[f"n_cells_{col[5:]}"] = ok.sum(axis=1)
        st[f"used_{col[5:]}"] = st.id.isin(df.columns)
    res = pd.DataFrame(out).sort_index()
    res.index.name = "date"; res = res.reset_index()
    res["date"] = res.date.dt.strftime("%Y-%m-%d")
    res.to_parquet(OUT, index=False); st.drop(columns="cell").to_csv(LST, index=False, lineterminator="\n")
    print(f"{len(res)} days -> {OUT.name}; stations used: " +
          ", ".join(f"{c}: {int(st[c].sum())}" for c in st.columns if c.startswith("used_")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
