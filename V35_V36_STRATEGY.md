# V35â€“V36: experimental-response retrieval and population views

This experiment follows V32â€“V34 without modifying their frozen code or results. V24 remains the published reference (overall 0.14299785521745675). The user requested both tests and a new upload-ready package. A package is not evidence that the private leaderboard score improves.

## Step 1 â€” inspect the existing pipeline

Reviewed source construction in `v30_methods.py`, paired-response features in `v32_methods.py`, V24 teacher/baseline logic, V16 integer count generation, the V34 held-family cases, and the original streaming writer and official packaging adapter. The main gaps addressed here are detailed experimental-response conditioning and supervision from independent populations. V34 lacked measured response tokens; V32/V33 used compact response programs. This experiment has no response-program decoder.

## Step 2 â€” use the existing measured response bank

Sources are full-target K562 and RPE1 pseudobulk plus H1, HCT116, HEK293T and three CD4 summary contexts. Pseudobulk is used as experimental response evidence, never represented as integer single cells. Missing genes remain absent. The present HepG2/Jurkat tables have zero rows; they are not claimed as training data. No new external source or private challenge response was obtained in this run.

V35 uses every gene measured in all allowed sources for its particular fitting fold. Output dimension can differ by fold. Other genes retain V24's transfer-effect values exactly; compositional normalization and integer sampling can still change their emitted counts. Control-rank coordinates exclude the 300 panel genes and use allowed controls only. Gene and perturbation features come from the previously verified public STRING embedding artifact.

## Step 3 â€” retrieve and correct V24

For each query perturbation, retrieve up to four exact/nearest STRING targets from each allowed source, capped at 32 tokens. Each token contains a measured gene-level response, source/control context coordinates, target features and candidate similarity. The initial prior accounts for family, source and candidate counts. A small learned attention network adjusts it based on recipient state and intervention. Each donor's directly targeted gene is excluded from that token; the query's own gene is also masked.

A shared gene decoder takes retrieved effect mean/spread/sign evidence, V24 effect, gene and target features, and control state. It predicts at most Â±0.5 log2FC correction. Its output layer starts at zero. Unsupported transfer-effect outputs are exact V24 fallback. The own target remains zero in the transferable effect array; the established emitter applies the explicit 90% knockdown.

Eight fixed training epochs, 160 steps per epoch, balanced biological-family sampling, AdamW, and target-identity 80/20 fitting/validation partitions are specified in the frozen plan. Epoch selection uses source-validation labels only. For every fitting recipient, both its complete family and the outer test family are excluded from V24 and retrieved response evidence. Query responses from other allowed contexts are legitimate conditioning information; target-held validation is not a claim that the target has never been experimentally observed anywhere.

## Step 4 â€” population-view training

V36 starts from its matching V35 fold. Public raw cells are available for H1, HCT116, HEK293T, CD4 D2 Rest and CD4 D2 Stim48. Only permitted fitting families are loaded. Each eligible source target produces eight pairs of disjoint treated populations (16â€“32 cells per view) and disjoint matched-lane control samples. Views can overlap across repetitions; each pair is disjoint. Sequencing is destructive: no individual cell is assigned an observed before/after counterpart.

One predicted mean is supervised against both independent population effects. A replicate-agreement weight reduces the influence of unstable gene changes. A second head learns bounded changes in the gene-wise coefficient of variation, with additional control variance/dropout features. Training balances summary-response and population-view losses. Recipient treatment mean or variance never appears in query conditioning. This is an independently implemented moment-summary adaptation, not the full P2P set encoder, an exact reproduction, or a guarantee that its paper's gains transfer.

Compare V36 means using the unchanged V24 emitter first. A separate arm changes only gene-wise empirical residual weights, while retaining the 0.65 Poisson component and dual mean constraints. This separates a mean improvement from a dispersion change. Null response/dispersion is zero and retains the original emitter; this approach does not claim to solve V24's pre-existing false-DE null artifact.

## Step 5 â€” test before packaging

1. Complete H1, complete Orion (both HCT116 and HEK293T), and complete CD4 exclusions.
2. Compare learned retrieval with V24, response identities shuffled within source, and fixed retrieval weights. These are inference ablations, not independently retrained shuffled models.
3. Independently replay all saved summary scores in float64; audit every source/family/target partition; independently implement the network forward in NumPy.
4. Evaluate both models on all 44 previously spent public raw mean cases.
5. Evaluate the same fixed 20 treatment count cases and five control nulls used for V34. D1 treatment cases contain 200 cells, so they are lower-power diagnostics; other populations contain 400.
6. Use the frozen release comparison rule. Explicitly report if strict improvement requirements fail. Do not open D3 or claim new blind validation. Keep V24 unchanged.
7. Fit the selected method on permitted full public sources, make 24 production preflight groups, then generate all 900 groups in official gene/context/target order.
8. Validate 360,000 integer-count cells Ã— 18,533 genes, 400 per target/context, finite nonnegative counts, row libraries, density and uniqueness. Invoke the unchanged official `vcc prep` after lossless bounded-memory staging, validate the wrapper, and record SHA-256. Do not upload automatically.

## Reproduce

Use Python3.12 and the pinned `requirements-v35-science.txt`; the exact executed environment is recorded in `experiments/v35_response_retrieval/analysis/runtime.json`. The competition environment is used for official packaging only. Frozen original output directories are immutable: use a new namespace for reruns after results exist. Dependencies and input hashes are recorded in the plan. Large arrays/models are ignored in Git; metadata, receipts, code and reports are committed.

```powershell
$py = 'C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe'
& $py -m unittest test_v35_v36 -v
& $py run_v35_v36.py freeze
& $py select_v35_v36_release.py freeze
& $py run_v35_v36.py train35
& $py run_v35_v36.py summary35
& $py audit_v35_v36.py 35
& $py run_v35_v36_cached.py verify
& $py run_v35_v36_cached.py train36
& $py run_v35_v36.py summary36
& $py audit_v35_v36.py 36
& $py audit_v35_v36_cache.py
& $py run_v35_v36.py mean
& $py run_v35_v36.py null
& $py run_v35_v36.py treatment
& $py audit_v35_v36_counts.py
& $py select_v35_v36_release_corrected.py select
# Fit V35 all first; V36 uses it if selected.
& $py run_v35_v36.py train35 --outer all
& $py run_v35_v36_cached.py train36 --outer all
& $py build_v35_v36_corrected_submission.py fit
& $py build_v35_v36_corrected_submission.py preflight
& $py build_v35_v36_corrected_submission.py generate
.\.competition_env\Scripts\python.exe package_v35_v36_submission.py
```

The corrected selector maps the frozen scorer's actual `common_expression_ratio` field to the existing expression-error comparison. The original selector is retained and a separate pre-outcome amendment records its correction. Numeric thresholds and scientific selection remain unchanged. Commands that freeze plans must not be run again in the original result directory.

## Primary research

- [PT-RAG](https://arxiv.org/abs/2603.07233): context-dependent differentiable retrieval over experimental perturbation evidence. We adapt the principle with STRING features and a small residual predictor; we do not use its STATE backbone/GenePT/Gumbel sampler or assert equivalent results.
- [P2P](https://arxiv.org/abs/2609.34391): population views share a repeatable perturbation effect while differing in sampling variation. Our two-view mean/dispersion supervision uses this motivation, with different architecture, sources and validation protocol.

Both papers' reported benchmark gains are separate from our VCC task. The public development results here cannot identify what private leading teams use or forecast their relative score.
For an outcome-known numerical rerun without overwriting any original plans or results:

```powershell
& $py reproduce_v35_v36.py
```

This creates a timestamped directory below `experiments/replays`. To resume one stage, pass the existing replay path with `--directory` and its `--stage`. The replay runner help and syntax were checked; a second full numerical training run is not claimed.

A separate normalization-cache amendment preserves the frozen original sampler. Before any V36 fit, both toy and actual raw-subset checks confirmed bit-identical means/variances. The cache stores only deterministic float64 per-cell normalization. Same cell rows, random-number order, labels, loss and epochs are retained. The initial pre-fit preparation log remains available; it contains no model score.

## Exploratory component extension

After both summary and44 raw-mean results and partial primary count outcomes were visible, one additional arm was frozen: V35 means with the fixed V36 dispersion head. It runs on the identical20cases, with unchanged control rows, seeds and variance weights. This tests whether variation helps when the stronger mean is retained. It is not a new blind or preregistered primary comparison. Original primary selection is preserved; a separate decision records whether fixed noninferiority/DE-improvement requirements permit adoption. No other hyperparameter search is introduced.

```powershell
& $py run_v35_v36_hybrid.py freeze
& $py run_v35_v36_hybrid.py run
& $py run_v35_v36_hybrid.py select
```

If this combination is selected, train V36 on all permitted sources for its variance head and retain the full-source V35 mean model. The artifact manifest records both checkpoints.

## Exploratory correction for discrimination drift

The full primary count tests showed lower PDS and higher corrected bulk-expression error despite improved CPM means. A retrospective all300-prediction diagnostic measured a target-shared correction component accounting for53–57% of correction energy in Orion contexts. Two fixed quarter-strength transforms were then frozen: target-center the V35 residual with exact support/own-gene and box constraints, or use the genuine-minus-shuffled retrieval difference. Both rerun all44 raw mean cases and the identical20countcases. This is post-outcome exploratory work, not new confirmation. Centering constrains log2-effect corrections across the given target panel, not raw counts or a universal biological law. Null response still uses the unchanged V24 emitter.

```powershell
& $py run_v35_v36_contrasts.py freeze
& $py run_v35_v36_contrasts.py run
& $py audit_v35_v36_contrasts.py
& $py select_v35_v36_contrasts_corrected.py
```

The conservative selection rule prefers count/mean noninferiority with some measured gain, and retains strict primary failure/confirmation flags. Original decisions are preserved. No further tuning grid, fresh reserve, upload or private-score claim is introduced.

The first contrast selector stopped on undefined HCT116 NMAE (fewer than10 reference DE genes). Its frozen original code and failure log are retained. The separate corrected selector applies the same original primary aggregation: omit undefined contexts within a study, require every study resolved, and never replace missing values by zero. No prediction, case, numeric threshold or scientific result changed. The full run uses one fixed training seed and one fixed emission seed per case; small gains are not statistical evidence of leaderboard improvement.

Production preflight rejected object-string label arrays before any cells were generated. The original generator, initial model arrays and metadata were preserved. `build_v35_v36_corrected_submission.py` encodes canonical CSV labels as NumPy Unicode and never unpickles the initial object arrays; both numeric effect and dispersion arrays were compared exactly after reencoding. Fresh production fits should use this adapter. An independent production audit verifies the actual V24 parent, official axes, finite bounds, own exclusion, exact fallback and genuinely changed predictions.

The360,000-row `prediction_control_rows.csv` stays local and ignored, consistent with previous submissions. Its exact sampler and seeds are committed; a final validation receipt records its SHA-256 and checks every original/predicted library. Group statistics, plans, model manifests and audits remain tracked.
