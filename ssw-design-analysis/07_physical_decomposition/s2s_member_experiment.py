#!/usr/bin/env python3
"""
s2s_member_experiment.py
========================
A NATURAL EXPERIMENT INSIDE OPERATIONAL ENSEMBLES: MEMBERS THAT DO AND DO NOT
REVERSE THE VORTEX FROM THE SAME INITIAL STATE.

Plan approved 2026-10-02 (user: "Please find a way to address the items 1-6!";
standing approval of new analyses). Committed BEFORE any member-level split by
reversal was computed. The per-member reforecast files on disk had been used only
for event-based tests (P', regional, lead-lag), never split by member reversal.

WHY
  The decisive causal evidence (SNAPSI) imposes two observed Northern Hemisphere
  SSWs. Every reforecast start of the ten S2S systems (4,500 starts, 1991-2024
  hindcast years, December-March) shares one initial state across its members;
  some members spontaneously reverse the 10 hPa, 60N wind and others do not. That
  is a matched comparison like SNAPSI's -- same initial conditions, SSW or no SSW
  -- over hundreds of winters, in ten other models, with SSWs that arise on their
  own rather than being imposed.

DATA (03_data_ingestion, per system: s2s_<c>_u10_60N[_long], _psl_cap[_short],
_t2m_regions[_short]; the P' confirmatory systems plus ECMWF; model versions as P')
  u(10 hPa, 60N) leads 1-34; polar-cap mean sea-level pressure leads 1-34;
  northern-Eurasian 2 m temperature leads 1-33.

MEMBERS
  Reversing member: westerly at leads 1-2, first westerly-to-easterly change of
  u(10 hPa, 60N) at a lead k_m in 3..20. Non-reversing: westerly at leads 1-2 and
  westerly throughout leads 3..20. Others (easterly at start) are dropped.
  Mixed start: at least two reversing and two non-reversing members.
  Anchor k* = median onset lead of the reversing members (rounded down).
OUTCOMES (window means, leads k*+8 .. k*+25 for the NAM proxy, k*+8 .. k*+24 for
  temperature; a start enters only if the window lies within the data, i.e.
  k* <= 9; secondary windows k*+8 .. k*+20 with k* <= 14)
  A: polar-cap NAM proxy = -(psl_cap - its mean over all starts of the system with
     the same start month-day and lead); T: northern-Eurasian temperature anomaly,
     same construction. Units: the system's s.d. of the window mean over all
     members of all starts (sigma_sys).
TESTS (pooled over systems and mixed starts; 2,000 bootstrap resamples of
  hindcast winters, clustering all starts and systems of a winter)
  N1 SHIFT: mean over starts of mean(reversing) - mean(non-reversing).
  N2 TRANSLATION: pooled within-start variance of reversing members about their
     own start mean over that of non-reversing members (each group's sum of
     squares over its degrees of freedom).
  N3 THRESHOLD: downward label (window mean < 0 and > 50% of days < 0) in each
     member; rate in each group; class contrast (DW minus NDW mean, pooled within
     group after removing start means) in each group; difference of contrasts.
  N4 STEP AT REVERSAL, INITIAL STATE FIXED: all members of all starts (westerly at
     leads 1-2); running variable = the member's minimum u over leads 3-20; key
     lead = lead of that minimum, <= 9; outcome window key+8 .. key+25; start fixed
     effects (outcome and running variable demeaned within start) ; dose slope per
     10 m/s and the coefficient of 1(u_min < 0) with separate slopes, as the
     observational test; covariate: the member's NAM proxy over leads 1-2.
  N5 BALANCE: NAM proxy over leads 1..k*-1 (before onset), reversing minus
     non-reversing -- how far members that will reverse already differ.
  N6 TEMPERATURE: N1-N3 for northern-Eurasian temperature.
  Descriptive: quantile shift (q 0.05..0.95) of reversing minus non-reversing,
  standardised within start; skewness of each group.
READING (fixed now): one shifted population predicts N1 < 0, N2 within about
  [0.9, 1.1] (no added spread), N3 contrast difference ~ 0, N4 dose > 0 with a step
  consistent with 0. Two kinds of response predict N2 > 1 (a mixture widens), a
  contrast that differs between the groups, or a step at reversal. A large N5
  means the comparison is not clean (members that reverse were already
  different) and is reported as such; N1 is then also reported net of N5.

Output: results/current/5_mechanism/s2s_member_experiment.json
"""
raise SystemExit("design registered 2026-10-02; implementation follows")
