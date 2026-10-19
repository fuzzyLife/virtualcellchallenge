# V69 live work log

User request: create a new submission after the POST_V68 literature review. Work is active; do not stop after research or a plan. No automatic upload or GitHub push. Commit locally as usual.

## Frozen science

- Root `v69_count_bridge.py`: independently implemented Gamma-Poisson posterior count coupling, beta/gamma dispersion-shape bridge, binomial thinning and Poisson additions. Exact zero-response count identity, varying final libraries, no global count-margin projection. Six numerical tests pass.
- `run_v69_experiment.py` frozen protocol committed as **383b827**; source fitting and production/audit helpers as **6fbea0f**. Plan SHA aac05c05ae0fabc7d0bddce03e7bbea59d69831823a048ce27038c77c8bac064. Do not modify frozen methods/runner/tests or old releases without preserving and disclosing corrections.
- Five source banks fitted: H1 25, HCT116234, HEK293T252, CD4 RestD2 34, Stim48D2 34 eligible target conditions. Whole-study shape exclusion. Baseline means are frozen public-fold V66 assay75; production exact published V66 arrays.
- Eligible candidates: mean_bridge, shape_bridge, hybrid25 (25% full shape bridge rows,75% unchanged V66 rows). First2 unique targets/context choose; remaining targets confirm. If none eligible, hybrid25 is predeclared **experimental fallback**. No confirmation-based retuning. Failed guards do not block the user-requested experimental artifact but must be disclosed.
- Five-arm labelled diagnostic splits each actual population into independent alternating fitting/scoring halves. Means are fitted, noisy diagnostics, not perfect oracles. Old-oracle has measured CPM+bulk, bridge-oracle uses bulk+shape, so not a pure shape comparison. Diagnostic never used for A/B/C fitting.

## Active work / runtime

Scientific Python: `C:/Users/animeshs/AppData/Local/Temp/vcc_v20_stack/env/Scripts/python.exe` with `-X utf8`. Official packaging uses `.competition_env/Scripts/python.exe`. UNC commands normally require `sandbox_permissions=require_escalated`; use login:false. Workspace has~1.7TB free, local temp~78GB free at check. V66 parent has2.368B nonzeros; memory-safe package staging is needed above2**31.

Diagnostic session24100 completed: raw-profile guard fails, oracle FID guard passes. Development **session16192** runs the two-worker context launcher, currently on the last CD4 stimulated context. Pilot tied-rank audit **session52309** may still be running; poll it once. Poll ongoing jobs with write_stdin. Source fit27036 finished; original null97137 finished; representative null32525 finished.

New `run_v69_cached.py` content-hashes complete CSR inputs and gene axis, caches identical DE/jackknife calculations with read-only arrays. Original numerical code remains unchanged. Cache self-test checks exact equality and input mutation keys. Use it, through `run_v69_parallel.py`, for development/confirmation to avoid repeated identical truth/reference computations. Parallel launcher uses two processes and per-context root experiment logs.

## Null sampling correction

Original frozen null used saved `null_template=x[left[:400]]`. HCT116/HEK293T real-count identities have578/1041 DE calls because prefix templates overrepresent batches. Original outcomes remain under analysis/null.

Separate **pre-additional-outcome** `representative_null_plan.json`/`run_v69_representative_null.py` freezes one depth-stratified sample of400 distinct controls over the full query-fit frame. All representative template/full bridge/near-zero bridge nulls have0 calls across5 contexts. V66 calls H1=1553,HCT1963,HEK1539,Rest1229,Stim2251; hybrid25=565,678,284,527,984 respectively. No seed tuning. This changes no selection rule. Pure null success is not a universal FPR guarantee.

## Matched control correction

Pre-development code review found original `shuffled_shape` is a full bridge population, so it is not a matched control for hybrid25. `run_v69_matched_control.py` and `matched_control_plan.json` add a quarter-mixture shuffled arm at the identical100/400 row positions. Frozen before any development outcome; no candidate or parameter is changed. Original code's beats_shuffle compares full bridge, so also report the **additional matched-control assessment**. The production model already labels hybrid25 unsupported because it retains null artifacts.

## Remaining execution

1. Wait for diagnostic session24100; inspect assessment. If missing metrics crash original assessor, preserve original failure and fix transparently, never substitute0. Expected HCT NMAE may be unavailable; do not invent it.
2. Run `run_v69_parallel.py development`. Then `run_v69_matched_control.py development`. Then original `run_v69_experiment.py select` (only once; refuses overwrite).
3. Run `run_v69_parallel.py confirmation`; then `run_v69_matched_control.py confirmation`; then `assess_v69_confirmation.py` (safe unavailable handling; original evaluator preserved). Original stricter metric guards may fail; disclose. Inspect matched control too.
4. Run `audit_v69_release.py results`. It verifies saved counts, independent LFC/BH/profile,8 manually reproduced rank coordinates/population, study traces/panel abstention. It also audits additional shuffled_hybrid rows automatically. Save source raw counts/arrays ignored.
5. `build_v69_submission.py fit`, `preflight`, `build --workers 2`. It loads exact parent V66 means/controls/provenance and V66 parent H5AD hash; quarter mixture reads original parent groups, other rows probabilistically transformed. Model status is experimental unless every applicable guard passes. **Review scripts before execution**; package path reference was fixed to v69_count_bridge.py.
6. Official `.competition_env/Scripts/python.exe package_v69_submission.py`. Lossless float32/int64 staging + unchanged official vcc prep. It needs local temp~30GB; do not run alongside the payload decompression audit.
7. `audit_v69_release.py payload`: complete payload/source comparison, all360k integer counts/libraries/identities, verifies all270k unchanged parent rows for hybrid, replays12 groups. Source count totals may fluctuate. No exact original-library equality assertion.
8. Run `document_v69_release.py` ONLY after complete audits and artifact. It writes V69_RESULTS.md and README. Matched-control and unavailable-metric reporting were added; helper not frozen yet. It also discloses the asymmetric oracle clipping/zero-activation constraints. It currently notes all confirmed failures and no hosted-score claim. Metric math is raw/equal-study, not anchor-scaled forecasts.
9. Add a release seal manifest/helper if useful; commit all new root codes/docs and compact JSON/CSV receipts, never ignored huge arrays/payloads or unrelated untracked files. Check prior scientific files remain unchanged. Give user ready path and Windows CMD submission command; label experimental if guards fail. V66 remains highest published0.150882188.

## Known code/interpretation considerations

- Original assessor directly compares aggregated NMAE, which might beNone if an entire study has too few real DE genes. Handle unavailable as unresolved/failed, not numeric0, if encountered.
- Frozen raw-profile selection is sensitive to increased sampling variance. Panel receipts ALSO retain signed jackknife noise-corrected expression numerators/denominators/rho. Document these descriptively; do not override selection after outcomes.
- Initial H1 bridge oracle improves FID with measured means, but raw profile error increases. In other contexts measured split-half means can be noisier than prediction. No broad conclusion from H1 alone.
- All actual400/200 row lists checked unique and treatment/template disjoint. No duplicate-cell manufacture. No new independent study or D3 opening.
- `package_v69_submission.py` reuses reviewed V66 lossless staging and official CLI; independent audit decompresses payload with `unpack_v24_official_payload.py`.
- `document_v69_release.py` is prepared but not run. `build_v69_submission.py` and `package_v69_submission.py` not yet run. Existing old payloads remain unchanged.
