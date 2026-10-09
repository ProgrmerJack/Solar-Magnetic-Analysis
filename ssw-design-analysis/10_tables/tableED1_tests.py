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
    M, N = J("snapsi_distribution_test.json"), J("snapsi_loeffel_test.json")
    WJ, AR = J("within_model_check.json"), J("snapsi_archetype_test.json")
    CS, OC = J("criterion_sweep.json"), J("obs_model_compatibility.json")
    FV, CM = J("forced_variance_ceiling.json"), J("class_model_comparison.json")
    HO, LL = J("s2s_heldout_test.json"), J("s2s_leadlag_test.json")
    HD = J("s2s_heldout_diagnosis.json")
    H26, VT = J("s2s_heldout2026_test.json"), J("vortex_threshold_continuity.json")
    SRF, RSD = J("shift_rule_forecast.json"), J("response_shape_dose.json")
    NX = J("s2s_member_experiment.json")
    RQ, TCV = J("snapsi_residual_quantiles.json"), J("criterion_transport_cv.json")
    SF, IC = J("state_forecast.json"), J("icon_event_heterogeneity.json")
    RCC = J("criterion_regional_contrast.json")
    GEN = J("snapsi_contrast_generality.json")
    PRO = J("prospective_next_ssw.json")
    DAV = J("davis2021_intervention_test.json")
    FRE = J("snapsi_free_natural_experiment.json")
    DAM = J("damip_class_law.json")
    DMM = J("damip_class_law_multimodel.json")
    CLG = J("class_law_general.json")
    MMT = J("mixture_moment_test.json")
    LCF = J("lens_confirmation.json")
    EOB = J("criterion_precip_contrast.json")
    RGM = J("snapsi_regime_precip.json")
    TCL = J("tc_class_law.json")
    EQV = J("equivalence_observational.json")
    ETR = J("early_transfer_retest.json")
    R2 = ("robustness (revision 2)", "committed before run", "4c8b3f2")
    mm = g(MMt, "multimodel", "confirmatory"); nr = g(SRt, "regions", "NEURASIA", "multimodel", "confirmatory")
    rows = [
        # (id, test, role, registration, commit, statistic, p)
        ("M-var", "variance ratio nudged/control, re-centred (SNAPSI NH)", "secondary", "in script, run same day", "0b6886d",
         f"{g(M,'pooled_all','variance_ratio')} {g(M,'pooled_all','variance_ratio_CI95')}", None),
        ("M-KS", "pure translation, Kolmogorov-Smirnov (SNAPSI NH)", "secondary", "in script, run same day", "0b6886d",
         "KS", g(M, "pooled_all", "pure_translation_KS_p")),
        ("N", "week-2 100 hPa coupling between members, nudged minus control (Fisher z)", "secondary",
         "in script with a design gate", "07324a4",
         f"dz {g(N,'primary','arm_contrast','delta_z_nudged_minus_control')}", g(N, "primary", "arm_contrast", "delta_p_two_sided")),
        ("J-P1", "within-model out-of-sample R2 added by the SSW, before onset (CMIP6)", "secondary", "in script, run same day", "0b6886d",
         f"{g(WJ,'P1 pre-onset','size_matched','event_specific')}", g(WJ, "P1 pre-onset", "size_matched", "p_null_ge_real")),
        ("J-P3", "within-model R2 added by the SSW, with post-onset stratosphere", "secondary", "in script, run same day", "0b6886d",
         f"{g(WJ,'P3 + post-onset stratosphere','size_matched','event_specific')}",
         g(WJ, "P3 + post-onset stratosphere", "size_matched", "p_null_ge_real")),
        ("O", "archetype pair: observed 2018-2019 difference within model pairings", "secondary", "in script, run same day", "79007c2",
         f"median percentile {g(AR,'pair_test','primary','median_obs_diff_pct')}; outside 95%: "
         f"{g(AR,'pair_test','primary','centres_outside_central_95')} of {g(AR,'pair_test','primary','n_centres')}", None),
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
        ("T2", "rank deficit, pre-onset starts that caught the SSW", "primary (revision)",
         "committed before data; matched null adopted after a review had emulated the polar-cap test (it raised p)", "57019f6, a201458",
         f"{g(PO,'outcomes','psl','T2_hits','mean_rank')}" if PO else "pending",
         g(PO, "outcomes", "psl", "T2_hits", "p_matched_null")),
        ("T2-miss", "rank deficit, pre-onset starts that missed the SSW", "secondary (revision)",
         "as T2", "57019f6, a201458",
         f"{g(PO,'outcomes','psl','T2_misses','mean_rank')}" if PO else "pending",
         g(PO, "outcomes", "psl", "T2_misses", "p_matched_null")),
        ("T2-d", "hit minus miss rank (polar cap)", "secondary (revision)", "as T2", "57019f6, a201458",
         f"{g(PO,'outcomes','psl','T2_hit_minus_miss','mean')} {g(PO,'outcomes','psl','T2_hit_minus_miss','ci95')}" if PO else "pending", None),
        ("G1", "split minus displaced surface response (ERA5)", "secondary (revision)", "committed before data", "57019f6",
         f"{g(VG,'tests','cls|nam_8_25','diff_split_minus_displaced')}" if VG else "pending",
         g(VG, "tests", "cls|nam_8_25", "p_perm_two_sided")),
        ("Q-P3", "CRPS skill with post-onset stratosphere: after SSWs vs ordinary days", "sensitivity", "post hoc", "this revision",
         f"{g(FR,'results','base | P3 + post-onset stratosphere','CRPSS_ssw')} vs "
         f"{g(FR,'results','base | P3 + post-onset stratosphere','null_ordinary_pooled','mean')}" if FR else "pending",
         g(FR, "results", "base | P3 + post-onset stratosphere", "null_ordinary_pooled", "p_null_ge_ssw")),
        # second revision (2026-09-30): robustness checks, each with its own registered
        # reading; not added to the primary Holm family
        ("T-a1", "criterion sweep, 54 surface versions: max abs z of SSW rate vs shifted null", *R2,
         f"max abs z {g(CS,'era5','S1_S2_conditions_1_2','max_abs_z')}" if CS else "pending",
         g(CS, "era5", "S1_S2_conditions_1_2", "p_max_abs_z")),
        ("T-a2", "criterion sweep, 54 surface versions: systematic excess (mean z)", *R2,
         f"mean z {g(CS,'era5','S1_S2_conditions_1_2','mean_z')}" if CS else "pending",
         g(CS, "era5", "S1_S2_conditions_1_2", "p_mean_z_excess")),
        ("T-a1-cal", "criterion sweep, surface versions: max abs z, calibrated reference", "robustness (revision 2)",
         "post hoc calibration of a registered test", "this revision",
         f"max abs z {g(CS,'era5','S1_S2_conditions_1_2','max_abs_z')}" if CS else "pending",
         g(CS, "era5", "calibration_posthoc", "c12", "p_max_abs_z_calibrated")),
        ("T-a2-cal", "criterion sweep, surface versions: mean z, calibrated reference (two-sided)", "robustness (revision 2)",
         "post hoc calibration of a registered test", "this revision",
         f"mean z {g(CS,'era5','S1_S2_conditions_1_2','mean_z')}" if CS else "pending",
         g(CS, "era5", "calibration_posthoc", "c12", "p_mean_z_calibrated_two_sided")),
        ("T-b", "observed statistics within 43-event CMIP6 draws (min-p combination)", *R2,
         f"min p {g(OC,'compatibility','pooled','min_p')}" if OC else "pending",
         g(OC, "compatibility", "pooled", "p_combined")),
        ("T-c", "identification audit: contrast / ceiling at SNAPSI-bounded damping (obs; CMIP6)", *R2,
         (f"{g(FV,'C4_identification_audit','datasets','observations','C4b_worst_case','ratio_contrast_over_ceiling')}; "
          f"{g(FV,'C4_identification_audit','datasets','cmip6','C4b_worst_case','ratio_contrast_over_ceiling')}") if FV else "pending",
         None),
        ("T-d", "two regimes beat continuous model after SSWs vs ordinary days (CRPS)", *R2,
         f"G {g(CM,'ssw','crps','G_M1_minus_M2')} vs {g(CM,'pseudo','crps','G_mean')}" if CM else "pending",
         g(CM, "pseudo", "crps", "p_pseudo_ge_ssw")),
        ("T-d-PH", "two regimes vs one skewed population, after SSWs vs ordinary days", "sensitivity", "post hoc",
         "this revision", f"G {g(CM,'posthoc','crps','PH1_G_M1s_minus_M2_ssw')}" if CM else "pending",
         g(CM, "posthoc", "crps", "PH1_p_pseudo_ge_ssw")),
        ("T-d-PH2", "two regimes vs continuous, calibrated on no-regime synthetic data", "sensitivity",
         "post hoc", "this revision", f"G {g(CM,'ssw','crps','G_M1_minus_M2')}" if CM else "pending",
         g(CM, "posthoc", "crps", "PH2_p_vs_negative_control")),
        ("T-d-PH4", "two-regime model skill over the shift: after SSWs vs ordinary days", "sensitivity",
         "post hoc", "this revision",
         (f"{g(CM,'posthoc','PH4_skill_vs_M0_crps','M2','ssw')} vs "
          f"{g(CM,'posthoc','PH4_skill_vs_M0_crps','M2','pseudo_mean')}") if CM else "pending",
         g(CM, "posthoc", "PH4_skill_vs_M0_crps", "M2", "p_pseudo_le_ssw")),
        ("T-d-PH3", "two-regime components after SSWs: mean separation (ordinary days)", "descriptive",
         "post hoc", "this revision",
         (f"{g(CM,'posthoc','PH3_ssw','separation')} ({g(CM,'posthoc','PH3_pseudo_mean_50sets','separation')})") if CM else "pending",
         None),
        ("HO1", "held-out 2023-24 SSWs: polar-cap rank deficit", "replication (revision 2)",
         "committed before data", "4c8b3f2",
         f"{g(HO,'HO1_polar_cap','mean_rank')} vs {g(HO,'HO1_polar_cap','null_mean')}" if g(HO, "HO1_polar_cap", "mean_rank") else "pending",
         g(HO, "HO1_polar_cap", "p")),
        ("HO2", "held-out 2023-24 SSWs: N-Eurasian temperature rank deficit", "replication (revision 2)",
         "committed before data", "4c8b3f2",
         f"{g(HO,'HO2_NEURASIA','mean_rank')} vs {g(HO,'HO2_NEURASIA','null_mean')}" if g(HO, "HO2_NEURASIA", "mean_rank") else "pending",
         g(HO, "HO2_NEURASIA", "p")),
        ("HD-2a", "ECMWF CY49R1 (model year 2025) on the 1998-2021 SSWs: rank deficit", "diagnostic",
         "post hoc (exploratory)", "this revision",
         f"{g(HD,'D2','ecmwf2025','mean_rank')} vs {g(HD,'D2','ecmwf2025','null_mean')}" if HD else "pending",
         g(HD, "D2", "ecmwf2025", "p")),
        ("HD-2b", "ECMWF CY49R1 minus CY47R3 rank, same 10 events", "diagnostic", "post hoc (exploratory)",
         "this revision",
         (f"{g(HD,'D2','ecmwf_same_events','diff_2025_minus_2022','mean')} "
          f"{g(HD,'D2','ecmwf_same_events','diff_2025_minus_2022','ci95')}") if HD else "pending", None),
        ("HD-2c", "CMA model year 2025 on the 1998-2021 SSWs: rank deficit", "diagnostic", "post hoc (exploratory)",
         "this revision",
         f"{g(HD,'D2','cma2025','mean_rank')} vs {g(HD,'D2','cma2025','null_mean')}" if HD else "pending",
         g(HD, "D2", "cma2025", "p")),
        ("HD-6", "held-out signs: likelihood ratio, forecasts as issued vs 1998-2021 bias", "diagnostic",
         "post hoc (exploratory)", "this revision",
         f"{g(HD,'D6','ratio_issued_over_recalibrated')}" if HD else "pending", None),
        ("L1", "100 hPa polar-cap height rank after SSWs", "diagnostic (revision 2)", "committed before data", "4c8b3f2",
         f"{g(LL,'tests','confirmatory','L1_z100_rank','mean_rank')} vs {g(LL,'tests','confirmatory','L1_z100_rank','null_mean')}" if LL else "pending",
         g(LL, "tests", "confirmatory", "L1_z100_rank", "p")),
        ("L2", "surface rank conditional on forecast 100 hPa anomaly", "diagnostic (revision 2)", "committed before data", "4c8b3f2",
         f"{g(LL,'tests','confirmatory','L2_conditional_surface_rank','mean_rank')} vs "
         f"{g(LL,'tests','confirmatory','L2_conditional_surface_rank','null_mean')}" if LL else "pending",
         g(LL, "tests", "confirmatory", "L2_conditional_surface_rank", "p")),
        ("HO2026", "held-out 4 March 2026 SSW, real-time forecasts: polar-cap rank deficit", "replication (revision 3)",
         "committed before data", "a511fbc",
         f"{g(H26,'HO2026_polar_cap','mean_rank')} vs {g(H26,'HO2026_polar_cap','null_mean')}" if H26 else "pending",
         g(H26, "HO2026_polar_cap", "p")),
        ("HO2026-T", "held-out 4 March 2026 SSW: N-Eurasian temperature rank deficit", "replication (revision 3)",
         "committed before data", "a511fbc",
         f"{g(H26,'HO2026_secondary_NEURASIA','mean_rank')} vs {g(H26,'HO2026_secondary_NEURASIA','null_mean')}" if H26 else "pending",
         g(H26, "HO2026_secondary_NEURASIA", "p")),
        ("HO-4", "four held-out SSWs (2023-2026) pooled: polar-cap rank deficit", "replication (revision 3)",
         "committed before data", "a511fbc",
         f"{g(H26,'four_events_polar_cap','mean_rank')} vs {g(H26,'four_events_polar_cap','null_mean')}" if H26 else "pending",
         g(H26, "four_events_polar_cap", "p")),
        ("U-dose", "AO days 8-52 per 10 m/s of minimum 10 hPa wind, 114 observed episodes", "revision-3 primary",
         "committed before outcomes", "2a0c40b",
         f"{g(VT,'observations','Y1_AO','E1_E2','dose_slope_per_10ms')} {g(VT,'observations','Y1_AO','E1_E2','dose_ci95')}" if VT else "pending",
         None),
        ("U-jump", "AO step at wind reversal, 114 observed episodes (Holm over 4 outcomes: 1.0)", "revision-3 primary",
         "committed before outcomes", "2a0c40b",
         f"{g(VT,'observations','Y1_AO','E1_E2','jump')} {g(VT,'observations','Y1_AO','E1_E2','jump_ci95')}" if VT else "pending",
         g(VT, "observations", "Y1_AO", "E1_E2", "jump_p_two_sided")),
        ("U-cmip6", "CMIP6 step at wind reversal (sigma), 5,726 episodes", "revision-3 secondary",
         "committed before outcomes", "2a0c40b",
         f"{g(VT,'cmip6','E1_E2','jump')} {g(VT,'cmip6','E1_E2','jump_ci95')}" if VT else "pending",
         g(VT, "cmip6", "E1_E2", "jump_p_two_sided")),
        ("V-0.10", "shift rule (SNAPSI) vs observed cold fortnights below 10th pct., 39 SSWs: binomial under rule", "revision-3 primary",
         "committed before run", "bf453bd",
         f"{g(SRF,'verification','V1','0.1','count')}/{g(SRF,'verification','V1','0.1','n')} vs p_rule {g(SRF,'verification','V1','0.1','p_rule')}" if SRF else "pending",
         g(SRF, "verification", "V1", "0.1", "binom_p_under_rule")),
        ("V-clim", "observed cold fortnights below 10th pct. vs climatology", "revision-3 primary",
         "committed before run", "bf453bd",
         f"{g(SRF,'verification','V1','0.1','count')}/{g(SRF,'verification','V1','0.1','n')} vs 0.10" if SRF else "pending",
         g(SRF, "verification", "V1", "0.1", "binom_p_under_climatology")),
        ("V-R1", "shift rule, strict independence (37 events, no SNAPSI events), q 0.10: binomial under rule", "revision-3 secondary",
         "committed before outcomes", "add7ba9",
         f"{g(SRF,'revision3','R1_strict_independence','0.1','count')}/{g(SRF,'revision3','R1_strict_independence','0.1','n')}" if SRF and SRF.get("revision3") else "pending",
         g(SRF, "revision3", "R1_strict_independence", "0.1", "binom_p_under_rule")),
        ("V-R2", "shift rule, 2.5th percentile, N Eurasia: binomial under rule", "revision-3 secondary",
         "committed before outcomes", "add7ba9",
         f"{g(SRF,'revision3','R2_R3_days8_24','NEURASIA','0.025','count')}/39 vs {g(SRF,'revision3','R2_R3_days8_24','NEURASIA','0.025','p_rule')}" if SRF and SRF.get("revision3") else "pending",
         g(SRF, "revision3", "R2_R3_days8_24", "NEURASIA", "0.025", "binom_p_under_rule")),
        ("V-R3", "shift rule, 2.5th percentile, high-latitude Europe: binomial under rule (rule rejected)", "revision-3 secondary",
         "committed before outcomes", "add7ba9",
         f"{g(SRF,'revision3','R2_R3_days8_24','HI_EUROPE','0.025','count')}/39 vs {g(SRF,'revision3','R2_R3_days8_24','HI_EUROPE','0.025','p_rule')}" if SRF and SRF.get("revision3") else "pending",
         g(SRF, "revision3", "R2_R3_days8_24", "HI_EUROPE", "0.025", "binom_p_under_rule")),
        ("V-R4", "shift rule, days 15-28, N Eurasia, q 0.10: binomial under rule", "revision-3 secondary",
         "committed before outcomes", "add7ba9",
         f"{g(SRF,'revision3','R4_weeks','days15_28','NEURASIA','0.1','count')}/39 vs {g(SRF,'revision3','R4_weeks','days15_28','NEURASIA','0.1','p_rule')}" if SRF and SRF.get("revision3") else "pending",
         g(SRF, "revision3", "R4_weeks", "days15_28", "NEURASIA", "0.1", "binom_p_under_rule")),
        ("V-R5", "second rule (CMIP6 x ERA5 event-free), N Eurasia q 0.10: binomial under rule", "revision-3 secondary",
         "committed before outcomes", "add7ba9",
         f"10/39 vs {g(SRF,'revision3','R5_second_predictor','p_rule2','0.1')}" if SRF and SRF.get("revision3") else "pending",
         g(SRF, "revision3", "R5_second_predictor", "verification", "0.1", "binom_p_under_rule")),
        ("V-R6", "shift rule OUT OF SAMPLE: 12 ERA5 SSWs 1941-58, N Eurasia q 0.10 (rule rejected)", "revision-3 primary",
         "committed before outcomes", "add7ba9",
         (f"{g(SRF,'revision3','R6_out_of_sample_1940_1958','primary','0.1','count')}/"
          f"{g(SRF,'revision3','R6_out_of_sample_1940_1958','primary','0.1','n')} vs 0.3243") if SRF and SRF.get("revision3", {}).get("R6_out_of_sample_1940_1958", {}).get("primary") else "pending",
         g(SRF, "revision3", "R6_out_of_sample_1940_1958", "primary", "0.1", "binom_p_under_rule")),
        ("U-ERA5", "AO step at reversal, ERA5 winds 1940-2025 (Holm over 4 outcomes)", "revision-3 secondary",
         "committed before data", "add7ba9",
         f"{g(VT,'era5_1940_2025','all','Y1_AO','E1_E2','jump')} {g(VT,'era5_1940_2025','all','Y1_AO','E1_E2','jump_ci95')}" if VT and VT.get("era5_1940_2025") else "pending",
         g(VT, "era5_1940_2025", "all", "E2_holm", "Y1_AO")),
        ("N-shift", "within-start ensemble comparison: shift, reversing minus non-reversing members (8-25 d, 40 starts)", "revision-3 primary",
         "committed before computation", "be804a7",
         f"{g(NX,'N1_shift_sigma','est')} {g(NX,'N1_shift_sigma','ci95')}" if NX else "pending", None),
        ("N-var", "within-start ensemble comparison: variance ratio (no widening)", "revision-3 primary",
         "committed before computation", "be804a7",
         f"{g(NX,'N2_variance_ratio','est')} {g(NX,'N2_variance_ratio','ci95')}" if NX else "pending", None),
        ("N-contrast", "within-start ensemble comparison: class contrast difference", "revision-3 primary",
         "committed before computation", "be804a7",
         f"{g(NX,'N3_threshold','contrast_diff')} {g(NX,'N3_threshold','contrast_diff_ci95')}" if NX else "pending", None),
        ("N-step", "within-ensemble step at reversal, initial state fixed (11,481 members)", "revision-3 primary",
         "committed before computation", "be804a7",
         f"{g(NX,'N4_step_initial_state_fixed','step_sigma')} {g(NX,'N4_step_initial_state_fixed','step_ci95')}" if NX else "pending",
         g(NX, "N4_step_initial_state_fixed", "step_p_two_sided")),
        ("RQ-A", "SNAPSI residual quantiles, polar-cap NAM, pooled: inside declared tolerance (+-0.25 sigma, +-0.05)?", "revision-3 primary",
         "committed before run", "113c7d3", g(RQ, "A_polar_cap_NAM_days8_25", "pooled", "reading") if RQ else "pending", None),
        ("RQ-T", "SNAPSI residual quantiles, N-Eurasian temperature, pooled", "revision-3 secondary",
         "committed before run", "113c7d3", g(RQ, "T_northern_Eurasia_days8_24", "pooled", "reading") if RQ else "pending", None),
        ("TCV", "shifted null, whole-winter cross-validation: observed vs expected downward count (retrospective)", "revision-3 primary",
         "committed before run", "bbfeb9c",
         f"{g(TCV,'primary_published_surface','observed_dw')} vs {g(TCV,'primary_published_surface','expected_dw')}" if TCV else "pending",
         g(TCV, "primary_published_surface", "p_calibration")),
        ("SF-F0", "state forecast vs climatology, 57 ERA5 SSWs, CRPS skill (retrospective CV)", "revision-3 primary",
         "committed before run", "4f5c53a",
         f"{g(SF,'E1_crps_skill','F3_vs_F0','all','skill')} {g(SF,'E1_crps_skill','F3_vs_F0','all','ci95')}" if SF else "pending", None),
        ("SF-F1", "state forecast vs fixed SNAPSI rule, CRPS skill", "revision-3 primary",
         "committed before run", "4f5c53a",
         f"{g(SF,'E1_crps_skill','F3_vs_F1','all','skill')} {g(SF,'E1_crps_skill','F3_vs_F1','all','ci95')}" if SF else "pending", None),
        ("SF-SSW", "SSW indicator beyond the continuous state, CRPS skill", "revision-3 primary",
         "committed before run", "4f5c53a",
         f"{g(SF,'E1_crps_skill','F4_vs_F3','all','skill')} {g(SF,'E1_crps_skill','F4_vs_F3','all','ci95')}" if SF else "pending", None),
        ("SF-ops", "calibrated operational ensembles vs state forecast, 17 events, CRPS skill", "revision-3 secondary",
         "committed before run", "4f5c53a",
         f"state vs calibrated {g(SF,'E2','F3_vs_cal','skill')} {g(SF,'E2','F3_vs_cal','ci95')}" if SF else "pending", None),
        ("ICON", "ICON event ensembles (18, each started at onset or the day before): between-event variance of mean responses, days 8-25 (bounds)",
         "revision-3 secondary", "committed before values read", "6ad6e8a",
         (f"{g(IC,'primary_registered','M2_between_event_variance','lower_bound_noise_sd_equals_daily_sd')}-"
          f"{g(IC,'primary_registered','M2_between_event_variance','upper_bound_noise_zero')}") if IC else "pending", None),
        ("RC", "N-Eurasian DW-NDW temperature contrast (days 8-24) vs matched shifted null, winter CV (power 0.30)",
         "revision-4 primary", "committed before run", "8c324ed",
         (f"obs {g(RCC,'primary_surface','cv','NEURASIA_d8_24','C_obs_K')} K vs null "
          f"{g(RCC,'primary_surface','cv','NEURASIA_d8_24','C_null_mean_K')} K") if RCC else "pending",
         g(RCC, 'primary_surface', 'cv', 'NEURASIA_d8_24', 'p_two_sided') if RCC else None),
        ("GEN-1", "back-shifted SSW-imposed contrast vs control, the 11 pairs with < 3 NDW members (tolerance 0.25 sigma); inconclusive",
         "revision-5 primary", "committed before run", "892377c",
         (f"{g(GEN,'G1_strongly_forced','strong_11','Delta_b','mean')} {g(GEN,'G1_strongly_forced','strong_11','Delta_b','ci95')}; "
          f"net of the spread change {g(GEN,'G1x_post_hoc_spread','strong_11','Delta_b_net_spread','mean')} "
          f"{g(GEN,'G1x_post_hoc_spread','strong_11','Delta_b_net_spread','ci95')} (post hoc); "
          f"deviations inflated by sqrt(1.5), about 0.4 of the control variance added: {g(GEN,'G1_strongly_forced','strong_11','Delta_b_k1.5','mean')} "
          f"{g(GEN,'G1_strongly_forced','strong_11','Delta_b_k1.5','ci95')} (not detected)") if GEN else "pending", None),
        ("GEN-2", "dose-response: back-shift difference on imposed shift, 36 pairs (slope per sigma)", "revision-5 secondary",
         "committed before run", "892377c",
         f"{g(GEN,'G2_dose_response','Delta_b_on_s_36','slope')} {g(GEN,'G2_dose_response','Delta_b_on_s_36','ci95')}" if GEN else "pending", None),
        ("GEN-3", "SH minor warming: nudged contrast vs location-scale null (and shift-only null)", "revision-5 primary",
         "committed before run", "892377c",
         (f"{g(GEN,'G3_third_event_SH','SH','resid_locscale','mean')} {g(GEN,'G3_third_event_SH','SH','resid_locscale','ci95')} "
          f"(shift only {g(GEN,'G3_third_event_SH','SH','resid_shift','mean')} {g(GEN,'G3_third_event_SH','SH','resid_shift','ci95')})") if GEN else "pending", None),
        ("GEN-4", "SNAPSI N-Eurasian DW-NDW temperature contrast vs matched shifted null (sigma; tolerance 0.25)", "revision-5 primary",
         "committed before run", "892377c",
         (f"{g(GEN,'G4_regional_temperature','resid','mean')} {g(GEN,'G4_regional_temperature','resid','ci95')}; "
          f"planted 0.5 sigma {g(GEN,'G4_regional_temperature','resid_plant_0.5sigma','mean')}") if GEN else "pending", None),
        ("GEN-5", "SNAPSI nudged-full (3-D stratosphere imposed) vs control, 9 pairs, 3 models: paired contrast difference; back-shift difference", "revision-5 primary",
         "committed before data", "892377c",
         (f"paired {g(GEN,'G5_nudged_full','control','paired_diff','mean')} (centre means ECMWF "
          f"{g(GEN,'G5x_post_hoc','control','centre_means','ECMWF','paired_diff')}, UKMO {g(GEN,'G5x_post_hoc','control','centre_means','UKMO','paired_diff')}); "
          f"back-shift {g(GEN,'G5_nudged_full','control','Delta_b','mean')} (centre means "
          f"{g(GEN,'G5x_post_hoc','control','centre_means','ECMWF','Delta_b')} to {g(GEN,'G5x_post_hoc','control','centre_means','Meteo-France','Delta_b')}); "
          f"2-3 centres, no interval; not decisive") if isinstance(g(GEN,'G5_nudged_full'), dict) else "pending", None),
        ("GEN-7", "why the zonal-pair elasticity is 1.14: slope with condition 1 only; calibration under the exact law", "revision-7 diagnostic",
         "committed before run", "a280a10",
         (f"condition 1 only {g(GEN,'G7_elasticity_diagnosis','a_condition1_only','slope')} {g(GEN,'G7_elasticity_diagnosis','a_condition1_only','ci95')}; "
          f"law calibration {g(GEN,'G7_elasticity_diagnosis','b_calibration_means_only','range95')}") if GEN else "pending", None),
        ("GEN-8", "shape of the outcome distribution, SSW imposed vs not, 36 SNAPSI pairs (AD permutation), with power against two-type mixtures", "revision-7 primary",
         "committed before run", "e78182a",
         (f"p {g(GEN,'G8_shape_with_power','AD_perm_p')}, max quantile difference {g(GEN,'G8_shape_with_power','max_abs_quantile_diff')}; "
          f"excluded: {len(g(GEN,'G8_shape_with_power','excluded_types') or [])} of 9 type cells (all six with separation 2-3, and 1.5 with half width)") if GEN else "pending",
         g(GEN, 'G8_shape_with_power', 'AD_perm_p') if GEN else None),
        ("GEN-8x", "test inversion of GEN-8 over two-Gaussian-type mixtures (separations 1-4, width ratios 0.5-2, any weights; common to the 36 pairs)", "revision-8 post hoc",
         "after review", "752fcb7",
         (f"{g(GEN,'G8x_post_hoc_inversion','n_rejected')} of {g(GEN,'G8x_post_hoc_inversion','n_family')} rejected; of the "
          f"{g(GEN,'G8x_post_hoc_inversion','signed','n_same_direction_half_shift')} able to produce half the imposed shift in its direction, "
          f"{g(GEN,'G8x_post_hoc_inversion','signed','n_same_direction_half_shift_rejected')} rejected; largest surviving implied shift "
          f"{g(GEN,'G8x_post_hoc_inversion','max_abs_implied_shift_among_survivors')} (observed {g(GEN,'G8x_post_hoc_inversion','observed_mean_abs_shift_sd')})") if GEN and g(GEN,'G8x_post_hoc_inversion','n_family') else "pending", None),
        ("GEN-9", "does the timing of the forced response explain the two-condition elasticity 1.14? time-varying-shift calibration under the exact law", "revision-8 primary",
         "committed before run", "752fcb7",
         (f"observed {g(GEN,'G9_timing_diagnosis','observed_slope_36')}; timing-only law {g(GEN,'G9_timing_diagnosis','c2','mean_slope')} "
          f"{g(GEN,'G9_timing_diagnosis','c2','range95')}; held out inside {g(GEN,'G9_timing_diagnosis','held_out_c2','share_inside')}; "
          f"{g(GEN,'G9_timing_diagnosis','reading')}") if GEN and g(GEN,'G9_timing_diagnosis','c2') else "pending", None),
        ("MM-1", "human forcing, three models (CanESM5, MIROC6, HadGEM3): central-European summer class law residual", "revision-7 primary",
         "committed before data", "a280a10",
         f"{g(DMM,'M1_CEUR_JJA_three_models','mean')} {g(DMM,'M1_CEUR_JJA_three_models','ci95')}; {g(DMM,'M1_CEUR_JJA_three_models','reading')}" if DMM else "pending", None),
        ("MM-2", "same: elasticity pooled over all model-region-years (spread barely changes: not discriminating)", "revision-7 primary",
         "committed before data", "a280a10",
         (f"{g(DMM,'M2_elasticity_all','slope')} {g(DMM,'M2_elasticity_all','ci95')}; expected under the law "
          f"{g(DMM,'post_hoc_calibration','expected_elasticity_under_law','range95')} (post hoc)") if DMM else "pending", None),
        ("MM-3", "same: winter cold classes, CanESM5 and MIROC6 (law residuals; registered reading)", "revision-7 secondary",
         "committed before data", "a280a10",
         ", ".join(f"{k.split('|')[0]} {k.split('|')[1]} {v['law']['mean']} ({v['reading']})" for k, v in (DMM or {}).get('per_model_target', {}).items() if 'DJF' in k and not k.startswith('HadGEM3')) if DMM else "pending", None),
        ("GL-1", "generality: law residual for temperature classes, human forcing, CanESM5 p1/p2, MIROC6, ACCESS-ESM1-5 (pooled, models equal)", "revision-8 primary",
         "committed before data", "80edfc5",
         (f"{g(CLG,'T1_by_variable','tas','pooled','mean')} {g(CLG,'T1_by_variable','tas','pooled','ci95')}; "
          f"share inside calibration {g(CLG,'T1_by_variable','tas','share_inside_calibration')} ({g(CLG,'T1_by_variable','tas','n_units')} cases); {g(CLG,'T1_by_variable','tas','reading')}") if CLG else "pending", None),
        ("GL-2", "same for precipitation classes (dry summers, wet winters), adding HadGEM3-GC31-LL", "revision-8 primary",
         "committed before data", "80edfc5",
         (f"{g(CLG,'T1_by_variable','pr','pooled','mean')} {g(CLG,'T1_by_variable','pr','pooled','ci95')}; "
          f"share inside calibration {g(CLG,'T1_by_variable','pr','share_inside_calibration')} ({g(CLG,'T1_by_variable','pr','n_units')} cases); {g(CLG,'T1_by_variable','pr','reading')}") if CLG else "pending", None),
        ("GL-3", "forced mixing of fixed response types: variance-mean curvature over natural/early-historical, historical, SSP5-8.5 (fixed types -1)", "revision-8 primary",
         "committed before data", "80edfc5",
         "; ".join(f"{v}: {g(CLG,'T3_by_variable',v,'n_mixing_excluded')} of {g(CLG,'T3_by_variable',v,'n_identified')} identified exclude mixing "
                   f"({g(CLG,'T3_by_variable',v,'n_units')} cases); {g(CLG,'T3_by_variable',v,'reading')}; calibrated against matched fixed-type mixtures "
                   f"(post hoc): {g(CLG,'T3_post_hoc_calibration_summary',v,'n_mixing_excluded_calibrated')} of {g(CLG,'T3_post_hoc_calibration_summary',v,'n_cases')} excluded; "
                   f"SUPERSEDED by MIX-1"
                   for v in ("tas", "pr")) if CLG else "pending", None),
        ("MIX-1", "direct moment test of forced fixed two-type mixing (validated before use; reanalysis of the GL-3 cases), Holm within variable", "revision-8 reanalysis",
         "method and gate fixed before model data", "70df9fa",
         "; ".join(f"{v}: {g(MMT,'summary',v+'|p_T','n_incompatible_holm')} of {g(MMT,'summary',v+'|p_T','n_cases')} incompatible "
                   f"({g(MMT,'summary',v+'|p_T','n_incompatible_raw')} at 5%)" for v in ("tas", "pr")) if MMT else "pending", None),
        ("CONF-1", "replication on untouched ensembles (CESM2 x2 forcing versions, MPI-GE; 1951-80, 2040-69, 2070-99): moment test of fixed two-type mixing, Holm, 36 cases per variable", "revision-8 primary",
         "committed before data", "d5096e4",
         "; ".join(f"{v}: {g(LCF,'summary','C1|'+v,'n_incompatible')} incompatible, {g(LCF,'summary','C1|'+v,'n_compatible')} compatible, "
                   f"{g(LCF,'summary','C1|'+v,'n_unresolved')} unresolved; {g(LCF,'summary','C1|'+v,'reading')}" for v in ("tas", "pr")) if LCF else "pending", None),
        ("CONF-2", "same: distribution (CDF) collinearity test", "revision-8 secondary", "committed before data", "d5096e4",
         "; ".join(f"{v}: {g(LCF,'summary','C2|'+v,'n_incompatible')} incompatible; {g(LCF,'summary','C2|'+v,'reading')}" for v in ("tas", "pr")) if LCF else "pending", None),
        ("CONF-3", "same: class law 2040-69 -> 2070-99 (pooled residual, models equal)", "revision-8 secondary", "committed before data", "d5096e4",
         "; ".join(f"{v}: {g(LCF,'summary','C3|'+v,'pooled','mean')} {g(LCF,'summary','C3|'+v,'pooled','ci95')}, {g(LCF,'summary','C3|'+v,'readings','holds')} of "
                   f"{g(LCF,'summary','C3|'+v,'n_cases')} hold; {g(LCF,'summary','C3|'+v,'reading')}" for v in ("tas", "pr")) if LCF else "pending", None),
        ("CONF-4", "same: class law 1951-80 -> 2040-69 (precipitation; temperature classes empty)", "revision-8 secondary", "committed before data", "d5096e4",
         (f"pr: {g(LCF,'summary','C4|pr','pooled','mean')} {g(LCF,'summary','C4|pr','pooled','ci95')}, {g(LCF,'summary','C4|pr','readings','holds')} of "
          f"{g(LCF,'summary','C4|pr','n_cases')} hold; {g(LCF,'summary','C4|pr','reading')}") if LCF else "pending", None),
        ("OBS-PR", "station rainfall (E-OBS) after the 39 observed SSWs: DW-NDW contrast vs matched null, Iberia days 1-25 (primary; Scandinavia, British Isles secondary)", "revision-8 primary",
         "committed before data", "9b0af12",
         (f"Iberia {g(EOB,'contrast','IBERIA_d1_25','C_obs')} vs null {g(EOB,'contrast','IBERIA_d1_25','C_null_mean')} mm/day, "
          f"p {g(EOB,'contrast','IBERIA_d1_25','p_two_sided')}, power {g(EOB,'contrast','IBERIA_d1_25','power_planted','detection_rate')}: "
          f"{g(EOB,'contrast','IBERIA_d1_25','reading')}; Scandinavia {g(EOB,'contrast','SCANDINAVIA_d1_25','C_obs')} vs "
          f"{g(EOB,'contrast','SCANDINAVIA_d1_25','C_null_mean')} (p {g(EOB,'contrast','SCANDINAVIA_d1_25','p_two_sided')})") if EOB else "pending",
         g(EOB, 'contrast', 'IBERIA_d1_25', 'p_two_sided') if EOB else None),
        ("EQ-1", "equivalence (TOST, +-0.25 event-free s.d.) and Bayes factor (half-normal, tau 0.5 s.d.) for the observational class-contrast tests: northern-Eurasian temperature, Iberian rainfall (primary); post hoc relative to those tests",
         "revision-8 post hoc", "bounds, priors, readings committed before computation", "7bce848",
         ("; ".join(f"{k.split('|')[1]}: TOST p {g(EQV,'results',k,'p_TOST')}, equivalent only at +-{g(EQV,'results',k,'smallest_bound_sd')} s.d., "
                    f"BF01 {g(EQV,'results',k,'BF01')}: {g(EQV,'results',k,'BF_reading')}" for k in ("T|NEURASIA_d8_24", "PR|IBERIA_d1_25"))) if EQV else "pending", None),
        ("R6-RT", "re-test of the failed 1941-1958 transfer with an independent detector (JRA-3Q) and outcome (GHCN-Daily stations); gates: detector >= 80% of 1958-80 catalogue onsets within 3 d, stations vs ERA5 r >= 0.90",
         "revision-8 primary", "committed before data", "3a2950b",
         (f"detector gate {g(ETR,'gate_D','matched_within_3d')}/{g(ETR,'gate_D','n_catalogue')}: FAILED, primary not run; station gate r {g(ETR,'gate_O','r_1959_2022')}: passed; "
          f"ERA5 events with station outcome {g(ETR,'arms','S1','0.1','count')}/{g(ETR,'arms','S1','0.1','n')} below q 0.10 (p under rule {g(ETR,'arms','S1','0.1','binom_p_under_rule')})") if ETR else "pending",
         g(ETR, 'arms', 'S1', '0.1', 'binom_p_under_rule') if ETR else None),
        ("REG-1", "SNAPSI: share of the imposed SSW's rainfall change reproduced by the frequency of four pressure-defined regimes (Iberia primary; Scandinavia), days 8-24, 36 pairs", "revision-8 primary",
         "committed before data (window amended before any outcome)", "6c99ebc, 44ff51a",
         ("; ".join(f"{b.title()}: forced {g(RGM,'k4','O1',b,'forced_change_mean_per_pair')} mm/day per pair, R {g(RGM,'k4','O1',b,'R')} {g(RGM,'k4','O1',b,'R_ci95')}; "
                    f"{g(RGM,'k4','O1',b,'reading')}" for b in ("IBERIA", "SCANDINAVIA"))) if RGM else "pending", None),
        ("TC-1", "tropical cyclones (IBTrACS 1982-2024, 3,740 storms; observational, SST levels confounded with basin): fixed two-type mixing across terciles of genesis SST, 24-h intensification and LMI (Holm over 2); ground truth from the observable RI label", "revision-8 primary",
         "committed before data (gate v2 before data)", "d3c6ef1",
         (f"rho {g(TCL,'C1_mixing','y_ri','rho')} {g(TCL,'C1_mixing','y_ri','rho_ci95')} and {g(TCL,'C1_mixing','y_lmi','rho')} {g(TCL,'C1_mixing','y_lmi','rho_ci95')}; "
          f"incompatible with fixed types in both (registered label: {g(TCL,'reading_pair')}; does not establish one population); ground truth {g(TCL,'C3_ground_truth','agreement_C1_lmi')}; law RI {g(TCL,'C2_law','RI','reading')}, major {g(TCL,'C2_law','MAJOR','reading')}") if TCL else "pending",
         g(TCL, 'C1_mixing', 'y_lmi', 'p_T_holm') if TCL else None),
        ("GL-4", "the law under SSP5-8.5 2071-2100 vs historical 1985-2014 (precipitation classes; temperature classes empty)", "revision-8 secondary",
         "committed before data", "80edfc5",
         (f"{g(CLG,'T4_summary','holds')} of {g(CLG,'T4_summary','n_units')} hold, {g(CLG,'T4_summary','inconclusive')} inconclusive, "
          f"{g(CLG,'T4_summary','fails')} fail") if CLG else "pending", None),
        ("DAV", "January 2021, CESM2(WACCM6) SSW vs scrambled-stratosphere ensembles (21 members), 4 Jan initialisation: back-shift contrast difference; not rejected, low power",
         "revision-5 primary", "committed before data", "738a422",
         f"{g(DAV,'pairs','04jan','Delta_b')} {g(DAV,'pairs','04jan','ci95','Delta_b')}" if DAV else "pending", None),
        ("DAV-pool", "the same, pooled over the 4 Jan, 1 Feb and 8 Feb initialisations", "revision-5 secondary",
         "committed before data", "738a422",
         f"{g(DAV,'R4_pooled_mean_of_pairs','Delta_b','mean')} {g(DAV,'R4_pooled_mean_of_pairs','Delta_b','ci95')}" if DAV else "pending", None),
        ("F2", "law on never-analysed data: SNAPSI free arm, members developing the SSW vs not (9 models); location-scale residual",
         "revision-6 primary", "committed before data", "1d93ff6",
         f"{g(FRE,'F2_resid_law','mean')} {g(FRE,'F2_resid_law','ci95')}; {g(FRE,'F2_reading')}" if FRE else "pending", None),
        ("F3", "same: elasticity of the contrast with respect to the spread (one population 1; two types >= 2 s.d. apart <= 0.60)",
         "revision-6 primary", "reading fixed before data (ea87721, corrected 51f89fe)", "1d93ff6",
         f"{g(FRE,'F3_Delta_b_vs_law','slope')} {g(FRE,'F3_Delta_b_vs_law','slope_ci95')}; {g(FRE,'F3_Delta_b_vs_law','reading')}" if FRE else "pending", None),
        ("F4", "same: shape of standardized groups (Anderson-Darling, 2,000 permutations)", "revision-6 secondary",
         "committed before data", "1d93ff6",
         f"p {g(FRE,'F4_shape','AD_perm_p')}, max quantile difference {g(FRE,'F4_shape','max_abs_quantile_diff')}" if FRE else "pending",
         g(FRE, 'F4_shape', 'AD_perm_p') if FRE else None),
        ("F5", "same: balance before onset (sigma)", "revision-6 secondary", "committed before data", "1d93ff6",
         f"{g(FRE,'F5_pre_onset_diff','mean')} {g(FRE,'F5_pre_onset_diff','ci95')}" if FRE else "pending", None),
        ("LAW-pool", "law residual pooled over five intervention sets (equal weight)", "revision-6 primary",
         "committed before the free-arm data; the other four sets' data had been analysed", "1d93ff6",
         f"{g(FRE,'L_pooled_law','location_scale','mean')} {g(FRE,'L_pooled_law','location_scale','ci95')}" if FRE else "pending", None),
        ("GEN-6", "regional cold-outcome classes, 24 regions x 36 SNAPSI pairs: law residual", "revision-6 primary",
         "committed before run", "1d93ff6",
         (f"{g(GEN,'G6_regional_classes','resid_law','mean')} {g(GEN,'G6_regional_classes','resid_law','ci95')}; "
          f"shift only {g(GEN,'G6_regional_classes','resid_shift','mean')}") if GEN else "pending", None),
        ("DAM", "second domain: human forcing, HadGEM3-GC31-LL historical vs hist-nat, central-European JJA upper-tercile class: law residual",
         "revision-6 primary", "committed before data", "60f3d8b",
         (f"{g(DAM,'targets','CEUR_JJA_tercile','resid_law','mean')} {g(DAM,'targets','CEUR_JJA_tercile','resid_law','ci95')}; "
          f"rate {g(DAM,'targets','CEUR_JJA_tercile','rate_nat','mean')} -> {g(DAM,'targets','CEUR_JJA_tercile','rate_hist','mean')}") if DAM else "pending", None),
        ("DAM-2", "same, western Russia and central North America JJA, central Eurasia DJF (tercile)", "revision-6 secondary",
         "committed before data", "60f3d8b",
         ", ".join(f"{t.split('_')[0]} {g(DAM,'targets',t,'resid_law','mean')}" for t in ("WRUS_JJA_tercile", "CNAM_JJA_tercile", "CEAS_DJF_tercile")) if DAM else "pending", None),
        ("DAM-3", "second domain: shape test (AD permutation) and elasticity, central-European JJA tercile", "revision-6 secondary",
         "committed before data", "60f3d8b",
         (f"shape p {g(DAM,'targets','CEUR_JJA_tercile','shape','AD_perm_p')}; elasticity "
          f"{g(DAM,'targets','CEUR_JJA_tercile','elasticity_slope','slope')} {g(DAM,'targets','CEUR_JJA_tercile','elasticity_slope','ci95')}") if DAM else "pending", None),
        ("DAM-4", "second domain: 10th-percentile classes (law residuals; the estimator is unbiased when the law holds exactly, post hoc)", "revision-6 secondary",
         "committed before data", "60f3d8b",
         ", ".join(f"{t.split('_')[0]} {g(DAM,'targets',t,'resid_law','mean')}" for t in ("CEUR_JJA_decile", "WRUS_JJA_decile", "CNAM_JJA_decile", "CEAS_DJF_decile")) if DAM else "pending", None),
        ("PRO", "prospective: shifted null on every NH SSW with onset Nov 2026 - Mar 2036 (CRPS, days 8-52 NAM)", "prospective primary",
         "predictions committed before any onset", "d2b54a8 (amended before any onset: 4cfba00, b2416e0)",
         (f"{g(PRO,'scores','cumulative','n_scored') or 0} events scored") if PRO else "pending", None),
        ("S-dose", "two regimes vs continuous with the realised 100 hPa dose (G, p vs ordinary days); inconclusive", "revision-3 primary",
         "committed before run", "7d8ed6f",
         (f"G {g(RSD,'S1_P2_plus_dose','crps','G_ssw')} vs {g(RSD,'S1_P2_plus_dose','crps','G_pseudo_mean')}; "
          f"power vs planted regimes {g(RSD,'posthoc_power','power_planted_1sd_regimes')} (post hoc): underpowered") if RSD else "pending",
         g(RSD, "S1_P2_plus_dose", "crps", "p_pseudo_ge_ssw")),
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
