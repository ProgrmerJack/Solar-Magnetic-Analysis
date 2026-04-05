"""
00_preregistration.py — Analysis Pre-Registration Protocol Lock
================================================================
Defines ALL analysis parameters BEFORE any data is examined.
This file serves as the pre-registration document for the study:

"Solar magnetic avalanches and terrestrial cryospheric instability:
 Evidence for a stratospheric coupling pathway"

ALL decisions below were made prior to running any statistical tests.
Git history of this file provides the audit trail.
"""
from dataclasses import dataclass, field
from datetime import date
import json
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# STUDY METADATA
# ═══════════════════════════════════════════════════════════════════════════════

STUDY_TITLE = (
    "Solar magnetic avalanches and terrestrial cryospheric instability: "
    "Evidence for a stratospheric coupling pathway"
)
PREREGISTRATION_DATE = "2026-04-05"
ANALYSIS_LOCKED = True  # Set True BEFORE running any analysis

# ═══════════════════════════════════════════════════════════════════════════════
# PRIMARY HYPOTHESIS
# ═══════════════════════════════════════════════════════════════════════════════

PRIMARY_HYPOTHESIS = (
    "Strong geomagnetic disturbance events (Kp >= 5 or Dst <= -50 nT) are "
    "followed by a statistically significant increase in Swiss natural "
    "avalanche activity within a 5-to-21-day lag window, mediated by "
    "stratospheric chemical and dynamical perturbations."
)

SECONDARY_HYPOTHESES = [
    "Geomagnetic disturbances produce measurable polar ozone depletion "
    "within 0-14 days (upstream manipulation check).",

    "Post-disturbance periods show weakened polar vortex (reduced 10 hPa "
    "zonal wind) within the 5-21 day window.",

    "A fast-pathway (24-72h) direct association between geomagnetic storms "
    "and avalanche activity exists, distinct from the stratospheric pathway.",

    "Dose-response: stronger disturbances (lower Dst, higher Kp) produce "
    "larger avalanche activity anomalies.",
]


# ═══════════════════════════════════════════════════════════════════════════════
# EXPOSURE DEFINITION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ExposureDefinition:
    """Geomagnetic disturbance event definition — locked before analysis."""

    # Primary thresholds
    kp_threshold: float = 5.0           # Kp_max >= 5 (G1+ storm)
    dst_threshold: float = -50.0        # Dst_min <= -50 nT (moderate storm)

    # Event definition uses OR logic: Kp >= threshold OR Dst <= threshold
    logic: str = "OR"

    # Declustering: minimum gap between independent events
    washout_days: int = 10

    # Temporal aggregation for daily exposure
    kp_agg: str = "max"     # Max 3-hourly Kp in calendar day
    dst_agg: str = "min"    # Min hourly Dst in calendar day (most negative)

    # Season restriction
    season_months: tuple = (11, 12, 1, 2, 3)   # NDJFM
    season_label: str = "NDJFM"

    # Sensitivity thresholds (tested in secondary analyses)
    kp_sensitivity: tuple = (4.0, 6.0, 7.0)
    dst_sensitivity: tuple = (-30.0, -75.0, -100.0)


@dataclass(frozen=True)
class FastPathwayDefinition:
    """Short-term (24-72h) CME/particle impact pathway."""
    lag_start_days: int = 1
    lag_end_days: int = 3
    label: str = "fast_pathway_1_3d"


@dataclass(frozen=True)
class StratosphericPathwayDefinition:
    """Medium-term (5-21d) stratospheric coupling pathway."""
    lag_start_days: int = 5
    lag_end_days: int = 21
    label: str = "strat_pathway_5_21d"


EXPOSURE = ExposureDefinition()
FAST_PATHWAY = FastPathwayDefinition()
STRAT_PATHWAY = StratosphericPathwayDefinition()


# ═══════════════════════════════════════════════════════════════════════════════
# PRIMARY OUTCOME
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class OutcomeDefinition:
    """Primary and secondary outcome measures — locked before analysis."""

    # Primary: daily natural avalanche activity index (SLF)
    primary_variable: str = "aai_all_natural"
    primary_description: str = "SLF daily natural avalanche activity index (all types)"

    # Secondary outcomes (tested with Bonferroni correction)
    secondary_variables: tuple = (
        "dry_natural_size_1234",    # Dry natural avalanches (all sizes)
        "wet_natural_size_1234",    # Wet natural avalanches
        "natural_size_234",         # Natural avalanches size >= 2
        "size_1234",                # All avalanches (any trigger)
        "max_size",                 # Maximum avalanche size that day
    )

    # Negative-control outcome (should show NO signal if mechanism is real)
    negative_control: str = "slf_accidents_daily_count"
    negative_control_description: str = (
        "SLF daily avalanche accident count (human-exposure-biased)"
    )


OUTCOME = OutcomeDefinition()


# ═══════════════════════════════════════════════════════════════════════════════
# STATISTICAL MODEL SPECIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ModelSpecification:
    """Primary regression model — locked before analysis."""

    # Model family
    family: str = "NegativeBinomial"  # for overdispersed count data
    link: str = "log"

    # Primary estimand: rate ratio = exp(β) for post-event window indicator
    estimand: str = "rate_ratio"

    # Fixed effects
    fixed_effects: tuple = ("winter_id", "month")

    # Confound controls (time-varying)
    confounders: tuple = (
        "nao_daily",            # North Atlantic Oscillation
        "qbo_u50",              # Quasi-Biennial Oscillation (50 hPa)
        "ncep_500hpa_nh",       # 500 hPa geopotential height (NH mean)
        "ncep_slp_nh",          # Sea-level pressure (NH mean)
        "day_of_season",        # Seasonality within winter
        "day_of_season_sq",     # Quadratic seasonality
    )

    # Standard errors
    se_method: str = "cluster_robust"
    cluster_variable: str = "winter_id"

    # Significance
    alpha: float = 0.05
    multiple_testing_correction: str = "Bonferroni"
    n_secondary_tests: int = 5   # Number of secondary outcomes


MODEL = ModelSpecification()


# ═══════════════════════════════════════════════════════════════════════════════
# CHEMISTRY MANIPULATION CHECK
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ChemistryCheckSpec:
    """Superposed epoch analysis for upstream chemistry validation."""

    # MLS variables to test
    variables: tuple = (
        "mls_ozone_1p0hpa", "mls_ozone_2p0hpa", "mls_ozone_4p6hpa",
        "mls_ozone_10p0hpa",
        "mls_temp_1p0hpa", "mls_temp_2p0hpa", "mls_temp_10p0hpa",
    )

    # Epoch window
    epoch_before_days: int = 15
    epoch_after_days: int = 30

    # Pre-specified sub-windows for aggregate test
    test_windows: tuple = ((0, 7), (8, 14), (15, 21))

    # Inference
    n_permutations: int = 2000
    permutation_method: str = "season_matched_sham"

    # Decision gate: if none of the chemistry variables show significant
    # response in the 0-14 day window, the causal chain cannot be established
    decision_gate: bool = True


CHEMISTRY_CHECK = ChemistryCheckSpec()


# ═══════════════════════════════════════════════════════════════════════════════
# DISTRIBUTED LAG MODEL (STRATOSPHERIC PROPAGATION)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class DistributedLagSpec:
    """Distributed lag model for stratospheric vortex response."""

    response_variables: tuple = ("ncep_u10_polar", "ncep_t10_polar")
    max_lag_days: int = 30
    lag_step: int = 1
    hac_bandwidth: str = "auto"   # Newey-West automatic bandwidth
    confidence_level: float = 0.95


DISTRIBUTED_LAG = DistributedLagSpec()


# ═══════════════════════════════════════════════════════════════════════════════
# FALSIFICATION TESTS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class FalsificationSpec:
    """Pre-specified falsification tests — all must pass."""

    # 7a. Sham event permutation
    n_sham_sets: int = 1000
    sham_kp_max: float = 3.0       # Quiet day threshold
    sham_dst_min: float = -20.0    # Quiet day threshold
    sham_ssw_buffer_days: int = 15
    sham_percentile_threshold: float = 95.0  # Real must beat 95% of shams

    # 7b. Summer null test
    summer_months: tuple = (5, 6, 7, 8, 9, 10)

    # 7c. Negative-control region
    control_region: str = "Norway"
    control_dataset: str = "norway_avalanche.parquet"

    # 7d. Negative-control outcome
    control_outcome: str = "slf_accidents"

    # 7e. Leave-one-winter-out cross-validation
    cv_method: str = "leave_one_winter_out"


FALSIFICATION = FalsificationSpec()


# ═══════════════════════════════════════════════════════════════════════════════
# STUDY PERIOD
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class StudyPeriod:
    """Temporal scope of the analysis."""

    # Primary analysis: MLS + SLF overlap
    primary_start: str = "2004-11-01"
    primary_end: str = "2019-05-31"
    primary_winters: int = 15  # 2004/05 through 2018/19

    # Extended analysis: Kp/Dst + SLF (no chemistry)
    extended_start: str = "1998-11-01"
    extended_end: str = "2019-05-31"
    extended_winters: int = 21

    # Winter definition
    winter_start_month: int = 11
    winter_end_month: int = 3


STUDY_PERIOD = StudyPeriod()


# ═══════════════════════════════════════════════════════════════════════════════
# DATA SOURCES
# ═══════════════════════════════════════════════════════════════════════════════

DATA_SOURCES = {
    "kp_index": {
        "file": "geomagnetic/kp_index.parquet",
        "role": "Primary exposure (3-hourly Kp, Ap)",
        "provider": "GFZ Potsdam / NOAA SWPC",
    },
    "dst_index": {
        "file": "geomagnetic/dst_index.parquet",
        "role": "Primary exposure (hourly Dst)",
        "provider": "Kyoto WDC",
    },
    "mls_ozone_polar": {
        "file": "atmospheric/mls_ozone_polar.parquet",
        "role": "Chemistry manipulation check",
        "provider": "NASA Aura/MLS v5",
    },
    "mls_temperature_polar": {
        "file": "atmospheric/mls_temperature_polar.parquet",
        "role": "Chemistry manipulation check",
        "provider": "NASA Aura/MLS v5",
    },
    "ncep_stratosphere": {
        "file": "atmospheric/ncep_stratosphere.parquet",
        "role": "Stratospheric dynamics (vortex)",
        "provider": "NCEP/NCAR Reanalysis",
    },
    "climate_indices": {
        "file": "atmospheric/climate_indices.parquet",
        "role": "Confound controls (NAO, QBO, ENSO)",
        "provider": "NOAA CPC / PSL",
    },
    "slf_activity": {
        "file": "cryosphere/slf_activity.parquet",
        "role": "PRIMARY ENDPOINT",
        "provider": "WSL/SLF Davos",
    },
    "slf_snow_events": {
        "file": "cryosphere/slf_snow_events.parquet",
        "role": "Individual event validation",
        "provider": "WSL/SLF Davos",
    },
    "slf_accidents": {
        "file": "cryosphere/slf_accidents.parquet",
        "role": "Negative-control outcome",
        "provider": "WSL/SLF Davos",
    },
    "flares": {
        "file": "solar/flares.parquet",
        "role": "Solar event context",
        "provider": "NASA DONKI",
    },
    "goes_xrs": {
        "file": "solar/goes_xrs.parquet",
        "role": "Solar X-ray context",
        "provider": "NOAA NCEI / GOES",
    },
    "omni_hourly": {
        "file": "solar/omni_hourly.parquet",
        "role": "Solar wind parameters",
        "provider": "NASA OMNIWeb",
    },
    "norway_avalanche": {
        "file": "cryosphere/norway_avalanche.parquet",
        "role": "Negative-control region",
        "provider": "NVE Varsom",
    },
    "ssw_catalog": {
        "file": "atmospheric/ssw_catalog.parquet",
        "role": "SSW event validation",
        "provider": "Butler et al. 2015",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════════

def export_preregistration(outpath: Path | None = None) -> dict:
    """Export all pre-registration parameters as JSON for reproducibility."""
    import dataclasses

    def dc_to_dict(obj):
        if dataclasses.is_dataclass(obj):
            return {k: dc_to_dict(v) for k, v in dataclasses.asdict(obj).items()}
        if isinstance(obj, tuple):
            return list(obj)
        return obj

    spec = {
        "study_title": STUDY_TITLE,
        "preregistration_date": PREREGISTRATION_DATE,
        "primary_hypothesis": PRIMARY_HYPOTHESIS,
        "secondary_hypotheses": SECONDARY_HYPOTHESES,
        "exposure": dc_to_dict(EXPOSURE),
        "fast_pathway": dc_to_dict(FAST_PATHWAY),
        "stratospheric_pathway": dc_to_dict(STRAT_PATHWAY),
        "outcome": dc_to_dict(OUTCOME),
        "model": dc_to_dict(MODEL),
        "chemistry_check": dc_to_dict(CHEMISTRY_CHECK),
        "distributed_lag": dc_to_dict(DISTRIBUTED_LAG),
        "falsification": dc_to_dict(FALSIFICATION),
        "study_period": dc_to_dict(STUDY_PERIOD),
        "data_sources": DATA_SOURCES,
    }

    if outpath is None:
        outpath = Path(__file__).resolve().parents[2] / "data" / "processed" / "preregistration.json"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(json.dumps(spec, indent=2, default=str), encoding="utf-8")
    print(f"Pre-registration exported to: {outpath}")
    return spec


if __name__ == "__main__":
    spec = export_preregistration()
    print(f"\n{'='*60}")
    print(f"STUDY: {STUDY_TITLE}")
    print(f"DATE:  {PREREGISTRATION_DATE}")
    print(f"LOCKED: {ANALYSIS_LOCKED}")
    print(f"{'='*60}")
    print(f"Exposure: Kp >= {EXPOSURE.kp_threshold} OR Dst <= {EXPOSURE.dst_threshold}")
    print(f"Washout: {EXPOSURE.washout_days} days")
    print(f"Season: {EXPOSURE.season_label}")
    print(f"Primary lag window: {STRAT_PATHWAY.lag_start_days}-{STRAT_PATHWAY.lag_end_days} days")
    print(f"Fast pathway window: {FAST_PATHWAY.lag_start_days}-{FAST_PATHWAY.lag_end_days} days")
    print(f"Primary outcome: {OUTCOME.primary_variable}")
    print(f"Model: {MODEL.family} ({MODEL.link} link)")
    print(f"SEs: {MODEL.se_method} by {MODEL.cluster_variable}")
    print(f"Alpha: {MODEL.alpha} ({MODEL.multiple_testing_correction})")
    print(f"Study period: {STUDY_PERIOD.primary_start} to {STUDY_PERIOD.primary_end}")
    print(f"Falsification tests: sham ({FALSIFICATION.n_sham_sets} sets), "
          f"summer null, {FALSIFICATION.control_region} control, "
          f"accident control, LOWO-CV")
    print(f"\n{len(DATA_SOURCES)} registered data sources")
