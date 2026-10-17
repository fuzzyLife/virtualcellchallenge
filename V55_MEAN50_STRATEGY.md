# V55: fixed V37 mean50 diagnostic submission

The user authorized reconsidering the previously unsubmitted V37 mean50 blend after reviewing correlations between public validation and official performance. Original V37 failed its REACH guard. Its original plan, selection and results remain unchanged. This is a separate, outcome-informed diagnostic; it is not an established better model.

The fixed response is float32(V24 + 0.5 * (submitted V36 retrieval_contrast25 - V24)). V24 residual10 emission, production template selection, own-target repression and seeds remain unchanged. No model retraining, coefficient search, oracle/private treatment values or new donor reserve is used.

The new plan is frozen before this run. Compare both arms on all44 available public populations at offsets0/550000,176 predicted-population comparisons in total; the40 original first-four/offset0 comparisons are reused by hash. These cover30 distinct context-target pairs, not176 independent biological experiments. CD4 Rest actual treatment populations contain200 cells; others400. Public fitting/validation biology and previous leaderboard outcomes have already been examined. D3 remains unopened.

Discrimination now ranks unique target identities within each context and seed. Repeated samples are averaged before ranking; legacy repeated-population PDS is retained separately. For DE summaries, first average repeated draws within context-target, then targets within context, contexts within study, and studies equally. Undefined NMAE remains absent. Expression ratios use the existing exact profile/jackknife scorer and retain warnings about tiny HCT116 denominators.

Diagnostic technical bounds are explicitly distinct from the original strict V37 release rule. Required: completed paired comparisons, legal counts/exact libraries; mean MSE/MAE<=1.002 baseline; PDS loss<=0.02; FID loss<=0.01; expression-ratio increase<=2%; available NMAE increase<=2%; unchanged null emitter. REACH and Jaccard outcomes must be disclosed. This authorization tests a trade-off and does not retrospectively convert the V37 failure into a pass. Expression/FID improvement in each seed and delete-target sensitivity are reported separately from technical eligibility. No coefficient changes or fallback arm selection follows these outcomes.

Independent checks replay saved-array fidelity, Jaccard and eligible NMAE, profile expression arithmetic, all30 delete-target expression comparisons and all5 exact null-emitter identities. Production checks independently reconstruct the float64 half-blend, official axes, exact fallback and target-own exclusion. Final delivery checks every library and row provenance, official packaging, and full decompressed CSR/axis equality. Large regenerable arrays, count files and packages stay local and outside Git.

Run scientific stages in the pinned Python3.12 environment; packaging uses the unchanged competition environment. Preserve original outputs; freeze refuses overwrite.

~~~bat
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" run_v55_mean50.py freeze
"%VCC_SCI_PY%" run_v55_mean50.py evaluate
"%VCC_SCI_PY%" run_v55_mean50.py assess
"%VCC_SCI_PY%" audit_v55_evaluation.py
"%VCC_SCI_PY%" run_v55_mean50.py fit
"%VCC_SCI_PY%" audit_v55_mean50.py production
"%VCC_SCI_PY%" run_v55_mean50.py preflight
"%VCC_SCI_PY%" run_v55_mean50.py generate
copy /b .\experiments\v55_mean50_diagnostic\analysis\selection.json .\experiments\v55_mean50_diagnostic\analysis\release_decision.json
.\.competition_env\Scripts\python.exe package_v55_mean50.py
"%VCC_SCI_PY%" audit_v55_mean50.py delivery
~~~

The compatibility receipt satisfies the frozen V50 packager's fixed release_decision filename. Its bytes equal the actual V55 selection receipt whose hash is bound in the model. The packaging adapter identifies inherited V50 template provenance fields and separately records the actual V55 generator/emitter/package adapter. Legacy generator cell IDs use a v36 prefix; these are unique synthetic identifiers and do not affect labels or counts. No upload is performed automatically. The finish_v55_mean50.py --wait helper chains assessment, independent audit, production/preflight/generation, byte-copy compatibility receipt, official packaging, delivery audit and documentation after evaluation completes. audit_v55_controls.py independently verifies the three original challenge-control SHA-256 values. Existing frozen destinations refuse overwrite; use a clean checkout/output namespace for a full numerical replay.
