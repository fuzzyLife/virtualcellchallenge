# V64 response-magnitude calibration: negative result

V64 separately tested the amplitude-only control discovered after V63 development. It multiplies V24 downstream LFC vectors by a teacher/V24 nonpanel vector-norm ratio clipped to 0..2, preserving signs and zeros. Own-target knockdown remains the existing separate 0.1 factor. This ratio is a model calibration feature, not measured knockdown efficacy. No neural training or new hidden-context intervention data were used.

The initial desired-mean diagnostic appeared to improve expression and PDS. Its 44 configurations included repeated target labels; the canonical unique-target reanalysis is saved separately. V64 was frozen as an explicitly outcome-informed follow-up, with unchanged public count controls and fixed scalar parameters.

## Integer-count result

The complete 60-comparison development screen fails. V24 versus magnitude, in raw public units:

| Metric | V24 | Magnitude |
| --- | ---: | ---: |
| PDS | 0.791667 | 0.763889 |
| Fidelity | 0.581703 | 0.580399 |
| Reach | 0.194389 | 0.229451 |
| Jaccard | 0.070027 | 0.070124 |
| NMAE error, lower is better | 0.912601 | 0.913525 |
| Expression ratio, lower is better | 37.471936 | 37.602475 |

Expression and discrimination fail, and the fixed allocation does not beat its scalar-shuffled control. The high expression-ratio magnitude is dominated by near-null HCT116 reference distance; it is not a private anchored score. V64 is rejected: no confirmation, production fit, count package or leaderboard submission is generated. V24 remains best published. This negative result also shows why an idealized mean-profile gain alone is insufficient to release a model.

## Scoring correction and provenance

The original frozen plan remains in initial_plan.json. While V64 development was running, before confirmation, a metadata audit found 44 predefined configurations span 30 unique target/context combinations and 34 distinct measured row sets, some overlapping. HCT116 and HEK293T reuse identical measured treated rows in two original configurations. Development already has unique labels and its verdict does not change. Any additional/confirmation panel must canonicalize by first case_index per literal target. The original outputs remain unchanged and corrected unique-panel results are separate.

Scripts: v64_response_magnitude.py, test_v64_magnitude.py, run_v64_counts.py, run_v64_development_parallel.py and audit_v64_unique_panels.py. Inputs, coefficients and scripts are hashed in experiments/v64_response_magnitude/analysis/plan.json; source folding and the count generator are inherited without change. Twenty identical V24 baseline statistics were copied with checksum receipts; these are computational reuse, not independent new measurements. All public biology was previously examined, and D3 remains unopened.

Use the existing scientific Python runtime:

```cmd
SCIENCE_PYTHON test_v64_magnitude.py
SCIENCE_PYTHON v64_response_magnitude.py freeze
SCIENCE_PYTHON amend_v64_scoring.py
SCIENCE_PYTHON v64_response_magnitude.py prepare
SCIENCE_PYTHON run_v64_development_parallel.py
SCIENCE_PYTHON audit_v64_unique_panels.py --version 64 --stage development
```

Freeze only in a new directory. No automatic upload or Git push.
