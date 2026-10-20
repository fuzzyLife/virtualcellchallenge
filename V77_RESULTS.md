# V77 results: common-response subtraction did not improve discrimination

The hypothesis was that nonspecific stress shared across many predicted interventions was obscuring target identity. V77 subtracted a leave-one-target-out common response at25%,50%,100%, or removed50% of its projection, preserving direct knockdown handling and bounded coordinate corrections. Source own-target coordinates were excluded from the common vector.

Before scoring,55 development target/context pairs were fixed: the first three legacy cases plus eight additional targets per context, chosen by a fixed hash from targets with at least100 cells and excluding every legacy case. The full655-pair truth bank was extracted, but the remaining600-pair confirmation was conditional on development success and was not run.

| Arm | Desired-profile error / V66 | PDS change |
|---|---:|---:|
| common25 |0.999670217| -0.003030303|
| common50 |0.999919642| -0.006060606|
| common100 |1.002159886| -0.007575758|
| projection50 |0.999470852| -0.003030303|

All four failed the prospective requirement of discrimination gain at least0.01 without a material profile regression. Expression-bin and panel-stratified gene-shuffled vectors were retained as negative controls. No remaining-target evaluation, count confirmation, production file or upload was generated.

This is reused public biology, not blind validation. The55-pair discrimination panel is more informative than the tiny legacy development panel; the profile means still use the first15 fixed cases. The evidence supports rejecting this particular subtraction, not the general claim that common biological responses never matter.

Code: `v77_common_response.py`, `prepare_v77_profiles.py`, `evaluate_v77.py`; protocol `V77_STRATEGY.md`; compact results `experiments/v77_common_response/analysis/development_results.json`.
