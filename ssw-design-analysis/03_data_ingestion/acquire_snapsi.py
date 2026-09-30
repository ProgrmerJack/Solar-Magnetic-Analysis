#!/usr/bin/env python3
"""
acquire_snapsi.py
=================
Acquires SNAPSI nudged-ensemble data — the experimental ground truth that lets
the pre-onset question be settled without relying on 43 observed events.

WHY SNAPSI
  Each SNAPSI case has three ensembles (GMD 15, 5073, 2022):
    free      atmosphere evolves freely
    nudged    stratosphere nudged to the OBSERVED evolution
    control   stratosphere nudged to CLIMATOLOGY
  nudged minus control is the stratospheric contribution BY EXPERIMENTAL DESIGN,
  with 50 members per centre per case instead of one realisation of history.
  That is the only way to resolve whether the pre-onset AO anomaly is a genuine
  precursor: observationally it is -0.82 on 43 events and every restricted
  subset's interval still contains it (see ../FINDING_design_robustness.md).

ARCHIVE LAYOUT — VERIFIED 2026-07-31, NOT GUESSED
  The previous version of this script probed "/badc/snapsi", got 404, and
  concluded SNAPSI was unreachable. That was a wrong path, not an outage. The
  real root is:

      https://data.ceda.ac.uk/badc/snap/data/post-cmip6/SNAPSI

  and the tree is:

      <centre>/<model>/<experiment>/s<YYYYMMDD>/<member>/<freq>/<var>/<grid>/<version>/<file>.nc

  Confirmed contents: 11 centres (CCCma, CNR-ISAC, ECCC, ECMWF, KMA,
  Meteo-France, NCAR, NOAA-GFDL, NRL, SNU, UKMO); experiments control/free/
  nudged; four start dates s20180125, s20180208, s20181213, s20190108; 50
  members; frequencies 6hr, 6hrPt, 6hrPtZ, 6hrZ.

  Variables actually present at 6hrPt (the useful frequency for surface work):
      hus mrso mrsos ps psl rlut siconc sithick snd snw ta tas tos ua uas va
      vas wap zg
  psl ~6.4 MB and tas ~6.8 MB per member-case; the 3-D fields (ua, zg, ta) are
  ~237 MB each and are deliberately NOT requested in bulk.

ACCESS
  Browsing and enumeration need NO credentials — the manifest below is built
  anonymously. Downloading a file redirects to auth.ceda.ac.uk/account/signin,
  so a **free CEDA account** is required for that step only
  (https://services.ceda.ac.uk/). This script does not and must not create one.

      set CEDA_USERNAME / CEDA_PASSWORD in the environment, or
      create ~/.netrc with a `machine dap.ceda.ac.uk` entry

  Without them the script still writes a complete manifest and downloads
  nothing, so acquisition becomes one command the moment credentials exist.

ENUMERATION STRATEGY
  Walking every member's variable/grid/version directory would be ~26,000
  requests. Instead the layout is DISCOVERED from the first member of each
  (centre, model, experiment, start-date) node, then the member id is
  substituted for the remaining members — and a random sample of the
  constructed URLs is VERIFIED against the archive. Any node whose sample fails
  verification is reported and its constructed rows are flagged
  `verified=False`, never silently trusted.

Outputs (this directory):
  snapsi_manifest.csv       every candidate file: URL, size, verified flag
  raw/snapsi/...            downloads, mirroring the archive layout (only with credentials)
"""
import hashlib
import os
import re
import sys
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "snapsi"
MANIFEST = HERE / "snapsi_manifest.csv"

BROWSE = "https://data.ceda.ac.uk"
DAP = "https://dap.ceda.ac.uk"
ROOT = "/badc/snap/data/post-cmip6/SNAPSI"

WANTED_FREQ = "6hrPt"
WANTED_VARS = ["psl", "tas"]        # 2-D only; 3-D fields are ~237 MB each
N_VERIFY = 3                        # constructed URLs sampled per node
TIMEOUT = 90


def session():
    s = requests.Session()
    s.headers["User-Agent"] = "ssw-design-analysis/1.0 (research)"
    u, p = os.environ.get("CEDA_USERNAME"), os.environ.get("CEDA_PASSWORD")
    if u and p:
        s.auth = (u, p)
        return s, True
    netrc = Path.home() / ".netrc"
    if netrc.exists() and "ceda.ac.uk" in netrc.read_text(errors="replace"):
        return s, True
    return s, False


def kids(s, path):
    """Immediate children of an archive directory (anonymous browse)."""
    try:
        r = s.get(BROWSE + path, timeout=TIMEOUT)
    except Exception:
        return []
    if not r.ok:
        return []
    out = []
    for h in re.findall(r'href="([^"]+)"', r.text):
        h = h.rstrip("/")
        if h.startswith(path + "/") and h.count("/") == path.count("/") + 1:
            n = h.split("/")[-1]
            if n not in out:
                out.append(n)
    return out


def leaf(s, path):
    """(filename, size_bytes) for the .nc in a version directory."""
    try:
        r = s.get(BROWSE + path, timeout=TIMEOUT)
    except Exception:
        return None, None
    if not r.ok:
        return None, None
    fn = re.findall(r"([A-Za-z0-9_\-]+\.nc)", r.text)
    if not fn:
        return None, None
    m = re.search(r"(\d+(?:\.\d+)?)\s*([KMG])i?B", r.text)
    size = None
    if m:
        size = float(m.group(1)) * {"K": 1e3, "M": 1e6, "G": 1e9}[m.group(2)]
    return fn[0], size


def main():
    s, creds = session()
    print(f"credentials present: {creds}"
          f"{'' if creds else '  -> manifest only, nothing will be downloaded'}")

    centres = [c for c in kids(s, ROOT) if c not in ("data", "post-cmip6", "SNAPSI")]
    if not centres:
        print(f"FAILED to list {ROOT} — archive layout may have changed")
        return 1
    print(f"centres: {centres}\n")

    rows, unverified_nodes = [], []
    for centre in centres:
        for model in kids(s, f"{ROOT}/{centre}"):
            for exp in kids(s, f"{ROOT}/{centre}/{model}"):
                for date in kids(s, f"{ROOT}/{centre}/{model}/{exp}"):
                    node = f"{ROOT}/{centre}/{model}/{exp}/{date}"
                    members = kids(s, node)
                    if not members:
                        continue
                    m0 = members[0]
                    for var in WANTED_VARS:
                        vp = f"{node}/{m0}/{WANTED_FREQ}/{var}"
                        grids = kids(s, vp)
                        if not grids:
                            continue
                        vers = kids(s, f"{vp}/{grids[0]}")
                        if not vers:
                            continue
                        vdir = f"{vp}/{grids[0]}/{vers[0]}"
                        fname, size = leaf(s, vdir)
                        if not fname:
                            continue
                        # substitute the member id into the discovered filename
                        made = []
                        for mem in members:
                            f2 = fname.replace(m0, mem)
                            p2 = (f"{node}/{mem}/{WANTED_FREQ}/{var}/"
                                  f"{grids[0]}/{vers[0]}/{f2}")
                            made.append((mem, p2, f2))
                        # verify a sample rather than trusting the pattern
                        import random
                        rnd = random.Random(20260731)
                        sample = rnd.sample(made, min(N_VERIFY, len(made)))
                        ok = 0
                        for mem, p2, f2 in sample:
                            d = "/".join(p2.split("/")[:-1])
                            got, _ = leaf(s, d)
                            if got == f2:
                                ok += 1
                        verified = ok == len(sample)
                        if not verified:
                            unverified_nodes.append(f"{centre}/{model}/{exp}/{date}/{var}")
                        for mem, p2, f2 in made:
                            rows.append({
                                "centre": centre, "model": model, "experiment": exp,
                                "start_date": date, "member": mem, "variable": var,
                                "frequency": WANTED_FREQ, "filename": f2,
                                "archive_path": p2,
                                "download_url": f"{DAP}{p2}?download=1",
                                "size_bytes_est": size, "verified": verified})
                    print(f"  {centre}/{model}/{exp}/{date}: {len(members)} members",
                          flush=True)

    if not rows:
        print("no files enumerated")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(MANIFEST, index=False, lineterminator="\n")

    tot = df["size_bytes_est"].sum() / 1e9
    print(f"\n{len(df):,} files enumerated across {df['centre'].nunique()} centres, "
          f"{df['experiment'].nunique()} experiments, {df['start_date'].nunique()} cases")
    print(f"  estimated total: {tot:.1f} GB for {WANTED_VARS} at {WANTED_FREQ}")
    print(f"  verified pattern: {int(df['verified'].sum()):,}/{len(df):,} rows")
    if unverified_nodes:
        print(f"  UNVERIFIED nodes ({len(unverified_nodes)}): {unverified_nodes[:5]}")
    print(f"  by experiment: {df.groupby('experiment').size().to_dict()}")
    print(f"-> {MANIFEST.name}")

    if not creds:
        print("\nNo CEDA credentials found — manifest written, nothing downloaded.")
        print("Register free at https://services.ceda.ac.uk/ then set")
        print("CEDA_USERNAME / CEDA_PASSWORD and re-run to fetch.")
        return 0

    RAW.mkdir(parents=True, exist_ok=True)
    got = 0
    for _, r in df[df["verified"]].iterrows():
        dest = RAW / r["archive_path"].lstrip("/")
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        resp = s.get(r["download_url"], timeout=600, stream=True, allow_redirects=True)
        if "auth.ceda.ac.uk" in resp.url:
            print("  authentication rejected — check credentials")
            return 1
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(1 << 20):
                fh.write(chunk)
        got += 1
        if got % 25 == 0:
            print(f"  downloaded {got}", flush=True)
    print(f"downloaded {got} files -> {RAW}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
