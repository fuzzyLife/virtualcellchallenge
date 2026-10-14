# V18: an independently pretrained response expert

Prepared 2026-10-13. This experiment tests a pretrained PIE expert against the released V16 response using complete emitted counts. The frozen comparison is `experiments/v18_pie_expert/analysis/plan.json` (SHA-256 `3cc7a58ac972c42c46f3969ec44faa05c729ba52552974c11a3793d427270bba`). Outcome: failed the frozen release gates; see `V18_RESULTS.md`. No V18 package was created.

## What is new

Use Arc's PIE K562-held-out checkpoint, trained on HepG2/Jurkat/RPE1. Its exact training metadata contains no overlap with our 300 targets. We use only the 6,290 competition genes measured in that training gene axis; V16 supplies all other coordinates. Gene descriptions alone do not establish trained output coverage.

Pinned author code: `ArcInstitute/pie@21a526add9ff3a81b4f962d9b16bbae570052c63`. Checkpoint: `arcinstitute/PIE_replogle_wdataset@055c7a2cabdaff5659121c64fa2b609c231fc5be/k562/best_auprc.ckpt`, SHA-256 `b9d5a0b3943666d892a56e7a9368f7850826b3f1fe6879e63a1a590f8add7d84`. See the upstream model/code licenses archived under `analysis/post_v17_research/`.

The model uses published ESM2 protein features, NCBI text tokens, STRING graph features, DepMap effects, perturbation descriptions and gene descriptions. Unknown context text is absent rather than a guessed cell identity. The 43,182,215,296-byte NCBI array is never downloaded in full: complete token sequences for 299 available panel targets occupy 744,243,200 payload bytes. Each bounded range has a receipt and the compact reconstruction is independently checked. We retain the published whole-array LFS hash as provenance, but do not claim to have computed that whole-array hash.

## Exact inference adaptation

All weights load strictly into the author's unchanged `PieModel`. We use Python 3.12, PyTorch 2.10 CUDA 12.6, full FP32, batch one, and 512-gene decoder chunks because the RTX 2070 SUPER lacks native BF16. The model's gene decoder has cross attention to shared latent memory, so narrowing the output axis does not add missing genes as zeros to its inputs.

Controls are means of **per-cell natural log1p(CP10k)**, using full source library totals before restricting genes. Taking log1p of mean CPM is a different operation and is not used. Benchmark inference reads only the left NTC controls; treatment counts enter only subsequent scoring.

All query datasets have novel names and every query target is absent from the complete training `pert_to_id`. Consequently the author's `serve_evidence` returns exactly zero donor and context evidence blocks. The adapter reproduces that restricted case directly, retaining the original training split hash and evidence key. A parity test uses the author's `serve_evidence`, `PieDataset` and `collate`; it would not justify skipping evidence for known training targets or contexts. This is an explicitly tested adapter, not a claim that the general official CLI needs no training assets.

Smoke test: VIM on public H1 controls, 6,290 outputs; finite predictions, strict weights and bitwise repeat passed. Two GPU forwards took 1.15 seconds with 466,024,448 peak allocated bytes. Raw predicted LFC ranged from -26.11 to 2.07; the predeclared clipping to [-3,3] and own-target exclusion are therefore important. No scoring outcome was used to choose this rule.

## Frozen comparison

Use the existing 14 development target/family pairs, two seeds (28 paired cases), with 400 distinct treated cells and 400 distinct query control templates per case. A separate right NTC population is the shared reference. Whole H1, HCT116 or HEK293T families are excluded from V16 response fitting in their respective benchmarks. Public biology was used in earlier work: this is algorithm verification, not blind independent biology.

Five arms use identical V16 `residual10` count generation and own-target knockdown:

1. Exact V16 comparator.
2. V16/PIE 75/25 LFC blend on supported coordinates.
3. Confidence blend with PIE weight `0.25 * clip((p_de - 0.5) / 0.5, 0, 1)`.
4. V16 with non-own L1 amplitude matched to the confidence blend.
5. Confidence blend using a frozen derangement of PIE target identities.

PIE `p_de` is not a BH q-value. Complete emitted counts are evaluated for raw PDS, corrected expression numerator/denominator, NMAE eligibility, fidelity, reach and Jaccard. All 300 panel genes are excluded from PDS; the perturbed gene is excluded from DE metrics. Raw errors are lower-is-better; these values are not the leaderboard's anchor-scaled scores.

The existing family-balanced and worst-family guards remain. At least three available metrics must improve, including a DE metric, and a finalist must beat both amplitude and shuffled controls on at least two of fidelity/reach/Jaccard. Missing HCT116 NMAE is reported as missing. Only a passing candidate may open the eight reserved target/family pairs (16 further paired cases), with the same once-only rule and no retuning. These are the corrected V17 verification identities; V17 never scored them.

V16's zero-effect count behavior remains unchanged. This response-only experiment does not claim to solve its null-DE or covariance limitations.

## Reproduction from Windows CMD

Run from the project root. Large matrices, weights, environments and regenerable NPZ files are deliberately ignored by Git; input hashes and compact records are retained.

```bat
.\.competition_env\Scripts\python.exe setup_v18_pie_runtime.py --install
.\.competition_env\Scripts\python.exe collect_v18_pie_assets.py
.\.competition_env\Scripts\python.exe audit_v18_assets.py
.\.competition_env\Scripts\python.exe test_v18_pie.py
.\.competition_env\Scripts\python.exe prepare_v18_experiment.py --inputs-only
"%TEMP%\vcc_v18_local\env\Scripts\python.exe" test_v18_pie.py --upstream
"%TEMP%\vcc_v18_local\env\Scripts\python.exe" infer_v18_pie.py smoke --family H1
.\.competition_env\Scripts\python.exe prepare_v18_experiment.py
"%TEMP%\vcc_v18_local\env\Scripts\python.exe" infer_v18_pie.py development --family H1
"%TEMP%\vcc_v18_local\env\Scripts\python.exe" infer_v18_pie.py development --family HCT116
"%TEMP%\vcc_v18_local\env\Scripts\python.exe" infer_v18_pie.py development --family HEK293T
.\.competition_env\Scripts\python.exe run_v18_experiment.py development --workers 4
.\.competition_env\Scripts\python.exe audit_v18_inputs.py
.\.competition_env\Scripts\python.exe audit_v18_experiment.py development
.\.competition_env\Scripts\python.exe assess_v18_experiment.py development
```

The freeze refuses to overwrite an existing plan. To reproduce a checked-out experiment, restore its exact scripts and plan, run `audit_v18_inputs.py --restore-missing` to regenerate ignored control inputs using the frozen NTC identities, and retain the pinned plan hash rather than refreezing against changed code. Inference and paired-case caches resume only with matching hashes; stage outputs are separate, so confirmation never overwrites development inference.

Primary sources: [PIE code](https://github.com/ArcInstitute/pie), [published checkpoint](https://huggingface.co/arcinstitute/PIE_replogle_wdataset), [PIE announcement](https://arcinstitute.org/news/pie). Broader alternatives and limitations are in `POST_V17_RESEARCH.md` and `POST_V17_NEXT_EXPERIMENT.md`.
