#!/usr/bin/env python3
"""
number_citations.py -- replace ^{key}^ citation keys in main.md with numbers in
Nature order (main text, then figure legends, then Methods) and write the
reference list. Every entry below was checked against Crossref or the publisher
page on 2026-09-25 (kretschmer18, huang21, kolstad20, garfinkel25 and dai25 on
2026-09-26). Run once the text is final; it rewrites main.md in place.
"""
import re
import sys
from pathlib import Path

REFS = {
    "baldwin21": "Baldwin, M. P. et al. Sudden stratospheric warmings. *Rev. Geophys.* **59**, e2020RG000708 (2021).",
    "sigmond13": "Sigmond, M., Scinocca, J. F., Kharin, V. V. & Shepherd, T. G. Enhanced seasonal forecast skill following stratospheric sudden warmings. *Nat. Geosci.* **6**, 98–102 (2013).",
    "karpechko17": "Karpechko, A. Y., Hitchcock, P., Peters, D. H. W. & Schneidereit, A. Predictability of downward propagation of major sudden stratospheric warmings. *Q. J. R. Meteorol. Soc.* **143**, 1459–1470 (2017).",
    "rao20": "Rao, J., Garfinkel, C. I. & White, I. P. Predicting the downward and surface influence of the February 2018 and January 2019 sudden stratospheric warming events in subseasonal to seasonal (S2S) models. *J. Geophys. Res. Atmos.* **125**, e2019JD031919 (2020).",
    "lu26": "Lu, R. & Rao, J. Sorting sudden stratospheric warmings with the downward tropospheric influence using ERA5 and CESM2-WACCM. *Atmos. Chem. Phys.* **26**, 3723–3742 (2026).",
    "nebel24": "Nebel, D. M. et al. The predictability of the downward versus non-downward propagation of sudden stratospheric warmings in S2S hindcasts. *Geophys. Res. Lett.* **51**, e2024GL110529 (2024).",
    "hitchcock14": "Hitchcock, P. & Simpson, I. R. The downward influence of stratospheric sudden warmings. *J. Atmos. Sci.* **71**, 3856–3876 (2014).",
    "hong26": "Hong, D.-C. et al. Surface impacts of Sudden Stratospheric Warmings (SSWs): comparison of 2018 and 2019 SSWs in SNAPSI experiments. Preprint at https://doi.org/10.5194/egusphere-2026-2798 (2026).",
    "leeRW25": "Lee, R. W., Charlton-Perez, A. J. & Lee, S. H. Stratospheric impacts on weather regimes following the 2018 and 2019 sudden stratospheric warmings. *Geophys. Res. Lett.* **52**, e2025GL115668 (2025).",
    "knight20": "Knight, J. et al. Predictability of European winters 2017/2018 and 2018/2019: contrasting influences from the Tropics and stratosphere. *Atmos. Sci. Lett.* **22**, e1009 (2021).",
    "white20": "White, I. P. et al. The generic nature of the tropospheric response to sudden stratospheric warmings. *J. Clim.* **33**, 5589–5610 (2020).",
    "loeffel26": "Loeffel, S. et al. Quantifying the tropospheric response to individual sudden stratospheric warmings revealed by an ensemble simulation strategy. *Weather Clim. Dynam.* **7**, 895–913 (2026).",
    "white19": "White, I. et al. The downward influence of sudden stratospheric warmings: association with tropospheric precursors. *J. Clim.* **32**, 85–108 (2019).",
    "bett23": "Bett, P. E., Scaife, A. A., Hardiman, S. C., Thornton, H. E. & Shen, X. Using large ensembles to quantify the impact of sudden stratospheric warmings and their precursors on the North Atlantic Oscillation. *Weather Clim. Dynam.* **4**, 213–228 (2023).",
    "hitchcock22": "Hitchcock, P. et al. Stratospheric Nudging And Predictable Surface Impacts (SNAPSI): a protocol for investigating the role of stratospheric polar vortex disturbances in subseasonal to seasonal forecasts. *Geosci. Model Dev.* **15**, 5073–5092 (2022).",
    "spaeth24": "Spaeth, J., Rupp, P., Garny, H. & Birner, T. Stratospheric impact on subseasonal forecast uncertainty in the northern extratropics. *Commun. Earth Environ.* **5**, 126 (2024).",
    "vitart17": "Vitart, F. et al. The Subseasonal to Seasonal (S2S) Prediction Project Database. *Bull. Am. Meteorol. Soc.* **98**, 163–173 (2017).",
    "feng25": "Feng, K. et al. Can stratospheric nudging improve surface predictability? Insights from the 2019 Southern Hemisphere sudden stratospheric warming. *npj Clim. Atmos. Sci.* **8**, 353 (2025).",
    "butler17": "Butler, A. H., Sjoberg, J. P., Seidel, D. J. & Rosenlof, K. H. A sudden stratospheric warming compendium. *Earth Syst. Sci. Data* **9**, 63–76 (2017).",
    "leeSH25": "Lee, S. H., Butler, A. H. & Manney, G. L. Two major sudden stratospheric warmings during winter 2023/2024. *Weather* **80**, 45–53 (2025).",
    "hersbach20": "Hersbach, H. et al. The ERA5 global reanalysis. *Q. J. R. Meteorol. Soc.* **146**, 1999–2049 (2020).",
    "rasp24": "Rasp, S. et al. WeatherBench 2: a benchmark for the next generation of data-driven global weather models. *J. Adv. Model. Earth Syst.* **16**, e2023MS004019 (2024).",
    "kretschmer18": "Kretschmer, M., Cohen, J., Matthias, V., Runge, J. & Coumou, D. The different stratospheric influence on cold-extremes in Eurasia and North America. *npj Clim. Atmos. Sci.* **1**, 44 (2018).",
    "huang21": "Huang, J., Hitchcock, P., Maycock, A. C., McKenna, C. M. & Tian, W. Northern hemisphere cold air outbreaks are more likely to be severe during weak polar vortex conditions. *Commun. Earth Environ.* **2**, 147 (2021).",
    "kolstad20": "Kolstad, E. W., Wulff, C. O., Domeisen, D. I. V. & Woollings, T. Tracing North Atlantic Oscillation forecast errors to stratospheric origins. *J. Clim.* **33**, 9145–9157 (2020).",
    "garfinkel25": "Garfinkel, C. I. et al. A process-based evaluation of biases in extratropical stratosphere–troposphere coupling in subseasonal forecast systems. *Weather Clim. Dynam.* **6**, 171–195 (2025).",
    "dai25": "Dai, Y., Hitchcock, P., Butler, A. H., Garfinkel, C. I. & Seviour, W. J. M. Assessing stratospheric contributions to subseasonal predictions of precipitation after the 2018 sudden stratospheric warming from the Stratospheric Nudging And Predictable Surface Impacts (SNAPSI) project. *Weather Clim. Dynam.* **6**, 841–862 (2025).",
    "charlton07": "Charlton, A. J. & Polvani, L. M. A new look at stratospheric sudden warmings. Part I: climatology and modeling benchmarks. *J. Clim.* **20**, 449–469 (2007).",
}
KEY = re.compile(r"\^\{([A-Za-z0-9,]+)\}\^")


def main(path):
    s = Path(path).read_text()
    if "[[REFERENCES]]" not in s:
        sys.exit("main.md has no [[REFERENCES]] placeholder: already numbered?")
    main_ = s[s.index("## Main"):s.index("## Methods")]
    legends = s[s.index("## Figure legends"):s.index("## Acknowledgements")]
    methods = s[s.index("## Methods"):s.index("## References")]
    order = []
    for part in (main_, legends, methods):
        for m in KEY.finditer(part):
            for k in m.group(1).split(","):
                if k not in order:
                    order.append(k)
    missing = [k for k in order if k not in REFS]
    if missing:
        sys.exit(f"no reference entry for {missing}")
    num = {k: i + 1 for i, k in enumerate(order)}

    def sub(m):
        ns = sorted(num[k] for k in m.group(1).split(","))
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(f"{ns[i]}–{ns[j]}" if j - i >= 2 else ",".join(map(str, ns[i:j + 1])))
            i = j + 1
        return "^" + ",".join(out) + "^"
    s = KEY.sub(sub, s)
    s = s.replace("[[REFERENCES]]", "\n".join(f"{num[k]}. {REFS[k]}" for k in order))
    Path(path).write_text(s, newline="\n")
    print(f"{len(order)} references numbered in order of first citation")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("main.md")))
