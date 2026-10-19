# V68 additional full-pipeline null audit

Code review after the initial H1 treatment results found a calibration gap in the fixed V68 protocol: its uniform-real null tests sampling of real cells, but does not exercise balancing against a zero-response target. The original results and methods are preserved. Before generating or viewing any null outcome, `primary_null_plan.json` freezes five extra complete-pipeline null comparisons.

For each context, optimize the same capped, depth-balanced weights against its frozen unperturbed bulk profile, with no own-gene edit. Sample 400 distinct controls and score against the same independent real reference, using the same null seed as the uniform-real arm. The candidate must additionally have no more than max(5,V63 null calls) significant genes in every context. This adds a requirement; it does not relax a failed guard or alter the candidate, regularization, solver or treatment predictions.

This addendum is disclosed as a code-review correction after some treatment outcomes, not a blind initial design. Its separate receipt binds the original plan and helper SHA-256. The tests still use previously studied public biology. No reserve data are opened.

```cmd
"%SCIENCE_PY%" -X utf8 run_v68_primary_null.py freeze
"%SCIENCE_PY%" -X utf8 run_v68_primary_null.py run
"%SCIENCE_PY%" -X utf8 run_v68_primary_null.py assess
```

The freeze refuses overwrite and requires both null-output folders to be empty. The run/check stages verify the pinned source before proceeding. These five added records are reported separately from the original ten null comparisons.
