#!/usr/bin/env python3
"""
tableED1_tests.py -- Extended Data Table 1: every test in the paper, whether it
was fixed before its data were analysed, where that is recorded, and its result.

Reads results/current/**.json (recomputes nothing except Holm's adjustment over
the family of registered primary tests that report a p value). Writes
09_figures/out/tableED1_tests.md and .csv.

Registration status, as the repository can show it:
  "committed before data"  -- the script with its design was committed before
                              the data it reads existed (git history)
  "committed before run"   -- committed before the confirmatory run, possibly
                              after the data were retrieved
  "in script, run same day"-- design written in the script docstring before the
                              run; first commit carries the result
  "post hoc"               -- added after the primary result (sensitivity)
"""
import csv
import glob
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "current"
OUT = ROOT / "ssw-design-analysis" / "09_figures" / "out"


def J(name):
    f = glob.glob(str(RES / "*" / name))
    return json.loads(Path(f[0]).read_text()) if f else None


def g(d, *path):
    for k in path:
        if d is None:
            return None
        d = d.get(k) if isinstance(d, dict) else None
    return d


def holm(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [None] * len(ps), 0.0
    for r, i in enumerate(idx):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = run
    return adj


def main():
    L, K = J("snapsi_selection_test.json"), J("snapsi_causal_effect.json")
    Q, P1 = J("forecast_value_test.json"), J("s2s_forecast_test.json")
    MMt, SRt = J("s2s_multimodel_test.json"), J("s2s_regional_test.json")
    SNr, ERr = J("snapsi_regional_test.json"), J("era5_regional_test.json")
    BL, RB = J("s2s_baseline_test.json"), J("s2s_multimodel_robustness.json")
    PO, VG = J("s2s_postonset_test.json"), J("vortex_geometry_test.json")
    CD, FR = J("snapsi_contrast_diagnosis.json"), J("forecast_value_robustness.json")
    mm = g(MMt, "multimodel", "confirmatory"); nr = g(SRt, "regions", "NEURASIA", "multimodel", "confirmatory")
    rows = [
        # (id, test, role, registration, commit, statistic, p)
        ("L", "paired DW-NDW contrast, nudged minus control (SNAPSI)", "primary", "in script, run same day", "0eda6d9",
         f"{g(L,'paired_NH','paired_difference')} {g(L,'paired_NH','paired_difference_CI95_centre_bootstrap')}", None),
        ("L-d", "nudged contrast vs empirical shifted-control null", "diagnostic", "post hoc", "this revision",
         f"{g(CD,'D2_resid_empirical','mean')} {g(CD,'D2_resid_empirical','ci95')}" if CD else "pending", None),
        ("R", "regional residual beyond the NAM, N Eurasia (SNAPSI)", "primary", "in script, run same day", "103a2c1",
         f"{g(SNr,'results','NEURASIA','R','mean')} {g(SNr,'results','NEURASIA','R','ci95')}", None),
        ("R-obs", "regional residual beyond the NAM, N Eurasia (ERA5)", "secondary", "in script, run same day", "103a2c1",
         f"{g(ERr,'regions','NEURASIA','R_mean')} K", g(ERr, "regions", "NEURASIA", "p_R_two_sided")),
        ("Q", "SSW-specific CRPS skill of event-aware forecasts (CMIP6, P2)", "primary", "in script, run same day", "0eda6d9",
         f"{g(Q,'tiers','P2 at-onset','ssw_specific_CRPSS')}", g(Q, "tiers", "P2 at-onset", "p_pseudo_ge_ssw")),
        ("P-D", "event discrimination, ECMWF (conditional null)", "discovery", "in script, run same day", "d9d8243",
         f"r {g(P1,'bins','short','discrimination','r_ensmean_vs_obs')}", g(P1, "bins", "short", "discrimination", "p_r_conditional")),
        ("P'-D", "event discrimination, nine-system mean", "secondary", "committed before run", "0eda6d9",
         f"r {g(mm,'D_r')}", g(mm, "D_p_conditional")),
        ("P'-H1", "outcomes low in ensembles, nine-system mean (polar cap)", "primary", "committed before run", "0eda6d9",
         f"rank {g(mm,'H1_mean_pit')} vs {g(mm,'H1_null_mean_pit')}", g(mm, "H1_p")),
        ("P'-H2", "shift-size correction improves CRPS", "primary", "committed before run", "0eda6d9",
         f"gain {g(mm,'H2_crps_gain')}", g(mm, "H2_p")),
        ("R-op", "N-Eurasian temperature low in ensembles (nine-system mean)", "primary", "committed before data", "103a2c1",
         f"rank {g(nr,'H1_mean_pit')} vs {g(nr,'H1_null_mean_pit')}", g(nr, "H1_p")),
        ("B1", "rank deficit vs all-winter baseline (polar cap)", "sensitivity", "post hoc", "this revision",
         f"{g(BL,'outcomes','psl','B0_mean_rank')} vs {g(BL,'outcomes','psl','B1_null_mean')}", g(BL, "outcomes", "psl", "B1_p")),
        ("B2", "rank deficit beyond vortex-state regression (polar cap)", "sensitivity", "post hoc", "this revision",
         f"residual {g(BL,'outcomes','psl','B2_residual_events')}", g(BL, "outcomes", "psl", "B2_p")),
        ("R3", "rank deficit, mixed model with date random effects (polar cap)", "sensitivity", "post hoc", "this revision",
         f"{g(RB,'outcomes','psl','R3_mixed_model','ssw_coef')} {g(RB,'outcomes','psl','R3_mixed_model','ssw_ci95')}",
         g(RB, "outcomes", "psl", "R3_mixed_model", "ssw_p")),
        ("T1", "rank deficit, starts after onset", "primary (revision)", "committed before data", "57019f6",
         f"{g(PO,'outcomes','psl','T1_post_onset','mean_rank')}" if PO else "pending",
         g(PO, "outcomes", "psl", "T1_post_onset", "p")),
        ("T2", "rank deficit, pre-onset starts that caught the SSW", "primary (revision)", "committed before data", "57019f6",
         f"{g(PO,'outcomes','psl','T2_hits','mean_rank')}" if PO else "pending",
         g(PO, "outcomes", "psl", "T2_hits", "p_matched_null")),
        ("G1", "split minus displaced surface response (ERA5)", "secondary (revision)", "committed before data", "57019f6",
         f"{g(VG,'tests','cls|nam_8_25','diff_split_minus_displaced')}" if VG else "pending",
         g(VG, "tests", "cls|nam_8_25", "p_perm_two_sided")),
        ("Q-P3", "SSW-specific CRPS skill with post-onset stratosphere", "sensitivity", "post hoc", "this revision",
         f"{g(FR,'results','base | P3 + post-onset stratosphere','CRPSS_ssw')}" if FR else "pending",
         g(FR, "results", "base | P3 + post-onset stratosphere", "null_ordinary_pooled", "p_null_ge_ssw")),
    ]
    fam = [i for i, r in enumerate(rows) if r[2].startswith("primary") and isinstance(r[6], (int, float))]
    adj = holm([rows[i][6] for i in fam])
    hp = {i: a for i, a in zip(fam, adj)}
    OUT.mkdir(parents=True, exist_ok=True)
    head = ["id", "test", "role", "registration", "commit", "statistic", "p", "p_Holm_primary"]
    with open(OUT / "tableED1_tests.csv", "w", newline="\n", encoding="utf8") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(head)
        for i, r in enumerate(rows):
            w.writerow([*r, round(hp[i], 4) if i in hp else ""])
    md = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for i, r in enumerate(rows):
        md.append("| " + " | ".join(str(x) if x is not None else "" for x in r) + f" | {round(hp[i], 4) if i in hp else ''} |")
    (OUT / "tableED1_tests.md").write_text("\n".join(md) + "\n", encoding="utf8", newline="\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
