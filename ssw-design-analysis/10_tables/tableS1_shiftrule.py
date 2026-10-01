#!/usr/bin/env python3
"""
tableS1_shiftrule.py -- Supplementary Tables 1-4: the shift rule in full.

Reads results/current/6_predictability/shift_rule_forecast.json (result V and its
revision-3 addendum). Recomputes nothing.

  S1  severity x region, days 8-24, 39 observed SSWs (SNAPSI probability, observed
      count/frequency with Wilson interval, binomial p under the rule and under
      climatology, risk ratios)
  S2  weeks 3-4 and 5-6
  S3  strict independence (without the SNAPSI events) and the second, CMIP6-based
      predictor
  S4  out of sample: SSWs of 1940-1958 in ERA5 (and from 1946; pooled with V1)
Output: 09_figures/out/tableS1_shiftrule.md
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "current" / "6_predictability" / "shift_rule_forecast.json"
OUT = ROOT / "ssw-design-analysis" / "09_figures" / "out" / "tableS1_shiftrule.md"
HEAD = "| region | window | q | rule p | observed | f [95%] | p (rule) | p (clim.) | RR obs | RR rule |"


def row(region, win, q, r):
    lo, hi = r["f_wilson95"]
    return (f"| {region} | {win} | {q} | {r['p_rule']} | {r['count']}/{r['n']} | {r['f']} [{lo}, {hi}] | "
            f"{r['binom_p_under_rule']} | {r['binom_p_under_climatology']} | {r['risk_ratio_observed']} | "
            f"{r['risk_ratio_rule']} |")


def block(title, d, region, win):
    out = []
    for q, r in d.items():
        if isinstance(r, dict) and "count" in r:
            out.append(row(region, win, q, r))
    return out


def main():
    d = json.loads(RES.read_text())
    r3 = d["revision3"]
    sep = "|" + "---|" * 10
    md = ["**Supplementary Table 1 | Severity and region, days 8–24 after 39 observed SSWs.**", "", HEAD, sep]
    for reg, v in r3["R2_R3_days8_24"].items():
        md += block("", v, reg, "8–24")
    md += ["", "**Supplementary Table 2 | Weeks 3–6.**", ""]
    for k, v in r3["R4_weeks"].items():
        sn = v["snapsi"]
        md.append(f"{k.replace('days', 'days ').replace('_', '–')}: SNAPSI {sn['n_pairs']} pairs, "
                  f"{sn['n_centres']} centres{'; ' + sn['status'] if sn.get('status') else ''}.")
    md += ["", HEAD, sep]
    for k, v in r3["R4_weeks"].items():
        for reg in ("NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"):
            if reg in v:
                md += block("", v[reg], reg, k.replace("days", "").replace("_", "–"))
    md += ["", "**Supplementary Table 3 | Strict independence and a second predictor (northern Eurasia, "
           "days 8–24).**", "", HEAD, sep]
    md += block("", r3["R1_strict_independence"], "NEURASIA (no SNAPSI events)", "8–24")
    md += block("", r3["R5_second_predictor"]["verification"], "NEURASIA (CMIP6 × ERA5 rule)", "8–24")
    reg5 = r3["R5_second_predictor"]["era5_free_regression"]
    md += ["", f"Second predictor: ERA5 event-free slope {reg5['slope_K_per_sd']} K per s.d. of the NAM "
           f"(residual s.d. {reg5['resid_sd']} K, {reg5['n_free']} windows); CMIP6 mean days 8–24 NAM after "
           f"{r3['R5_second_predictor']['cmip6_n_events']} SSWs {r3['R5_second_predictor']['cmip6_mean_nam_days8_24']}."]
    r6 = r3["R6_out_of_sample_1940_1958"]
    if r6.get("status"):
        md += ["", f"**Supplementary Table 4 | Out of sample, 1940–1958:** {r6['status']}."]
        OUT.write_text("\n".join(md) + "\n", encoding="utf8", newline="\n")
        print(f"-> {OUT.relative_to(ROOT)} (R6 pending)")
        return
    md += ["", "**Supplementary Table 4 | Out of sample: ERA5 SSWs of 1940–1958 (outcomes not examined "
           "before registration).**", "",
           f"Onsets: {', '.join(r6['era5_onsets_1940_1958'])}. Detector against the catalogue: "
           f"{r6['detector_vs_catalogue']['matched_within_3d']} of {r6['detector_vs_catalogue']['n_catalogue']} "
           "catalogued onsets matched within 3 days.", "", HEAD, sep]
    md += block("", r6["primary"], "NEURASIA 1940–58", "8–24")
    md += block("", r6["from_1946"], "NEURASIA 1946–58", "8–24")
    md += ["", "| pooled with V1 | q | observed | f [95%] | p (rule) | p (clim.) |", "|---|---|---|---|---|---|"]
    for q, r in r6["pooled_with_V1"].items():
        md.append(f"| NEURASIA | {q} | {r['count']}/{r['n']} | {r['f']} [{r['f_wilson95'][0]}, {r['f_wilson95'][1]}] | "
                  f"{r['binom_p_under_rule']} | {r['binom_p_under_climatology']} |")
    s = r6["secondary_negative_nam"]
    md += ["", f"Negative polar-cap NAM proxy, days 8–52: {s['count']}/{s['n']} against the CMIP6 0.74 "
           f"(p {s['binom_p_under_rule']}) and the period's climatological {s['climatological_base_rate']} "
           f"(p {s['binom_p_under_climatology']})."]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(md) + "\n", encoding="utf8", newline="\n")
    print(f"-> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
