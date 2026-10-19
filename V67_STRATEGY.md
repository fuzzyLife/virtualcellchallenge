# V67: intervention-network rewiring

V66 improved overall by only 0.000459 over V63, with reach offsetting five regressions. Its small local count test got five component directions wrong. This experiment changes the response hypothesis and uses the broader cohort first, rather than selecting another amplitude or blend on the same 30 pairs.

## New hypothesis

The 300 target genes are both interventions and measured outputs. Their cross-responses can define a feedback operator. Transfer may depend on how much the intermediate target genes are available in the recipient context, not just which source response has the nearest control profile.

For each allowed source family, take knockdown-adjusted target-to-target log2 responses, soft-threshold at two source standard errors, and construct S = I + response-transpose / log2(0.1). A regularized inverse and diagonal normalization provide an approximate interaction operator. Bound its spectral norm at 0.45. Solve source and query feedback resolvents, normalizing each intervention column by its diagonal to represent clamping its own gene. The query multiplies outgoing mediator edges by a fixed bounded saturation ratio derived from query and source basal CPM. Transport the source's centered downstream responses through the change in resolvent; subtract the unmodified source response, precision-pool supported changes and apply half of a correction bounded at one log2 unit onto V63.

This is a predictive approximation inspired by [modular response analysis](https://doi.org/10.1006/jtbi.2002.3096) and [CellOracle's signal propagation](https://www.nature.com/articles/s41586-022-05688-9). It is not an implementation of either paper, an identified causal graph, or a claim that 90% CRISPRi is an infinitesimal perturbation. The saturation law is a fixed hypothesis, not a recovered biological law. Perturbations of the 300-gene subnetwork can have unobserved mediators. Strong intervention, noisy feedback inversion and non-transcriptional targets can invalidate this model.

This differs from V25's covariance-response geometry, V30's STRING-to-response ridge model and V35's embedding-based response retrieval: the edges here come from experimentally measured cross-target effects, and inference re-solves a feedback system. It introduces no new neural network or pretrained weight. It uses existing verified source evidence; no new independent dataset is claimed.

## Controls and invariants

Only the fixed `network` arm is eligible. Compare unchanged V63 and V66, mediator-availability permutation, and an amplitude-only projection of the network correction with matching output support. Controls cannot be promoted to a submission after observing outcomes. All 300 panel output genes retain V63; they serve as mediators rather than being directly changed by the correction. Existing own-target suppression remains separate and occurs once during eventual cell generation.

Seven numerical tests check exact unchanged-context identity, zero-network identity, recovery of a known feedback system, mediator-chain propagation, bounded stability, missing-pair abstention including measured zeros, and poisoning of invalid response coordinates. An independent audit poisons all recipient-study treatment responses, uncertainties, validity masks, knockdown summaries and source profiles and requires every actual source-model field to remain bitwise identical. HCT116 and HEK293T are jointly excluded as Orion; all CD4 conditions are excluded together.

Unsupported measurements abstain. Downstream network outputs require valid measurements for every active intervention in that source. No output is imputed as a true zero response. Neutral target nodes with no measured intervention response remain in the identity operator and cannot add a fitted source response.

## Evaluation frozen before outcomes

Predict all 300 responses in every context before loading recipient treatment matrices for scoring. Use all existing 655 target/context pairs with at least 30 genuine cells across H1, HCT116, HEK293T, resting CD4 donor 2 and stimulated CD4 donor 2. These represent three studies. Compute intended bulk-profile error, CPM fold-change error and target discrimination, excluding own genes and excluding all panel genes from discrimination.

Whole-target SHA256 partitions are shared across contexts, and development/confirmation summaries are reported separately. The public biological datasets have been used before; these partitions are diagnostic, not newly blind data. A canonical-small target subset is also scored on the same full treatment populations, isolating target-set composition within this evaluation rather than reproducing the older sampled-cell test.

Require at least 1% lower equal-study profile error than BOTH V63 and V66; no study worse by more than 2%; PDS loss at most 0.005 versus either; lower profile error in both target partitions; and a better profile error than both negative controls. Only a passing primary proceeds to paired integer-cell tests with the unchanged residual10 emitter. All count gates, source hashes and exact formulas are in the immutable `analysis/plan.json`. No post-outcome parameter tuning. D3 remains unopened. A failed primary is documented and does not become a renamed submission.

## Reproduce

Use the recorded scientific Python environment, not the official packaging environment:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v67_network.py
"%SCIENCE_PY%" -X utf8 v67_intervention_network.py freeze
"%SCIENCE_PY%" -X utf8 audit_v67_boundaries.py
"%SCIENCE_PY%" -X utf8 v67_intervention_network.py prepare
"%SCIENCE_PY%" -X utf8 run_v67_experiment.py means
```

Freeze refuses overwrite. Completed per-context predictions and results resume only with matching hashes. Scripts remain in the root; compact receipts and reports are versioned; large arrays and runtime logs are ignored. No upload or Git push occurs automatically. Results and any stopping decision are recorded separately.
