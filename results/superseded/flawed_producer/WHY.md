# Results whose producing script is known to be wrong

Moved out of `results/current/` on 2026-09-17.

## `gate3_final.json`

Written by `06_simulation_validation/finalise_gate3.py`, which the project's own
`FINDING_gate3_contamination.md` line 11 labels **"(flawed)"**: its `excluded`
pseudo-onset draw leaves real-event-influenced days inside the baseline, so the
null centres near **+0.29 σ** instead of zero. The file records that directly —
`specs.harm6.excluded.null_mean = 0.2933`, `passes_strict = false`.

It was dangerous where it sat. Its name implies it is the final word on Gate 3,
it lived beside the correct artifact in `results/current/3_calibration/`, and the
two disagree. Anyone resolving the conflict by filename would have taken the
wrong one.

**The correct artifact is `results/current/3_calibration/gate3_clean_null.json`**,
produced by `gate3_clean_null.py` at N_BOOT=1000 / N_COVER=300 on the 43-event
catalogue: 3 harmonics null mean **+0.0052 [−0.0062, +0.0166]**, coverage
**0.930**, FPR **0.070**, `passes = true`. Gate 3 closes on that file.

The flawed script itself is still in the repo. `FAILURES.md` (2026-08-04) records
the rule this produced: when two scripts answer one question, delete or rename the
dead one — a file named `*_final.py` is not evidence that it is.
