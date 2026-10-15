"""Reproduce actual V26 calibration outcomes, README and artifact decision."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v26_methods import ROOT,A,FAMILIES,KINDS,load_plan,write
import build_submission as b

def main():
 plan=load_plan();audit=json.loads((A/'independent_audit.json').read_text());assert audit['passed'] and audit['plan_sha256']==b.sha256(A/'plan.json');stopfile=A/'calibration_stop.json';assert stopfile.exists(),'This reporter handles the calibration-rejected path; passing models need the locked reserve and full production audit';stop=json.loads(stopfile.read_text());key=stop['failed_family'];result=json.loads((A/f'null_{key}_result.json').read_text());assert not result['passed'];assert not (A/'development_evaluation.json').exists() and not (A/'reserve_execution_started.json').exists();models=[json.loads(p.read_text()) for p in sorted((A/'models').glob('*.json'))];tables=[pd.read_csv(A/f'null_{f}_scores.csv') for f in FAMILIES if (A/f'null_{f}_scores.csv').exists()];table=pd.concat(tables,ignore_index=True);assert audit['null_cases_checked']==table.family.nunique() and audit['decoder_models_checked']==len(models);rows=[]
 for _,r in table.iterrows():rows.append('|'+str(r.family)+'|'+str(r.candidate)+'|'+str(int(r.null_DE_calls))+'|'+str(int(r.real_null_DE_calls))+'|'+f'{r.detection_MAE:.6f}|{r.log_variance_MAE:.6f}|{r.mean_reference_L1:.8f}'+'|')
 training=json.loads((A/'training'/f'{key}.json').read_text());fit_n=sum(sum(r['fit_flags']) for r in training['original_inputs']);valid_n=training['rows']-fit_n;champion=json.loads((ROOT/'experiments/v24_cd4_response/analysis/leaderboard_result.json').read_text());write(A/'release_decision.json',{'released':False,'uploaded':False,'new_package_generated':False,'reason':'Conditional count decoder failed prospectively frozen real-control null calibration','failed_family':key,'failed_checks':stop['failed_checks'],'treatment_metrics_opened':False,'reserve_opened':False,'retained_version':24,'retained_entry_id':champion['entry_id'],'retained_overall':champion['overall'],'plan_sha256':b.sha256(A/'plan.json'),'calibration_stop_sha256':b.sha256(stopfile),'independent_audit_sha256':b.sha256(A/'independent_audit.json'),'assessment_recovery_sha256':b.sha256(A/'assessment_recovery.json'),'latent_diagnostic_sha256':b.sha256(A/'H1_latent_diagnostic.json'),'null_gene_diagnostic_sha256':b.sha256(A/'H1_null_gene_diagnostic.json'),'documenter_sha256':b.sha256(Path(__file__))});report=f'''# V26 results: neural count decoder

**No V26 package is released.** The conditional neural decoder fails the frozen {key} real-control calibration gate. Training and count generation completed for the executed pilot, but treatment outcomes and the donor reserve remain unopened. V24 stays the published reference at{champion['overall']}. This is a negative result for this source-trained decoder configuration, not a test of every neural count model.

## What actually ran

Plan SHA256:`{b.sha256(A/'plan.json')}`. Fitted {len(models)} prespecified neural models, each for12 epochs, using{training['rows']:,} sampled source cells ({fit_n:,} fitting/{valid_n:,} validation) across{len(training['training_contexts'])} training contexts, with{training['trained_output_genes']:,} learned output coordinates. The18533-gene head has explicit masks; unsupported coordinates retain empirical query variation. The complete held-out family is excluded from training, state feature selection, the response basis and the V24 mean predictor. Source data and original cell counts were independently checked.

The neural likelihood models16 latent dimensions and eight state dimensions, with a64/128-unit decoder and16 intervention descriptors. Model epochs are selected by held-cell validation ELBO, which is not a competition metric. The independent, conditional and training-state-shuffled arms share samples and architecture. All fitted weights, source identities, training-only group means and epoch traces are retained.

Executed {len(table)} paired count populations across{table.family.nunique()} frozen null cases. Each has400 actual controls,400 independent reference controls and400 distinct template controls drawn from the separate inference pool. Six population arms are evaluated. Target-descriptor shuffling is identical to conditional for a zero intervention and reuses that count population; this null screen does not establish target-specific treatment accuracy.

## Measured raw null diagnostics

DE calls use the existing competition-style CPU Wilcoxon/BH scorer. Detection and log-variance errors compare each generated population with the same real reference controls. These are raw diagnostic values, not anchor-scaled leaderboard scores.

|Context|Arm|Generated DE calls|Actual null DE calls|Detection MAE|Log-variance MAE|Mean CPM probability L1 vs V24|
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Failed checks: **{', '.join(stop['failed_checks'])}**. The mean contract, DE-call envelope, detection envelope and variance envelope are separate checks; see the exact Boolean outcomes in null_{key}_result.json. Returning V24 exactly is not treated as calibrated null performance.

The V26 conditional pipeline reduces null DE calls from821 to352 (57.1%) and improves log-variance error from0.8942 to0.6198. Detection error slightly worsens from0.06080 to0.06220; independent real-control detection error is0.01599. Three calibration gates fail; the mean contract passes at about0.000000952 L1 drift. These are one-case raw diagnostics, not treatment gains or a leaderboard prediction. The state-shuffled decoder has337 calls, so the call-count reduction alone does not establish that correct state conditioning is responsible.

The empirical/resampling comparator produces6764 calls and detection error0.15999. Its shared sampling/repair pipeline substantially changes observed-cell distributions. This warns against treating re-sampled observed count probabilities as a neutral empirical reference: the input counts already contain measurement noise. It does not isolate the effects of sampling, stratum constraints and swaps. New generator improvements are collective pipeline effects until separate ablations identify them.

The first numerical sampler designs were rejected on synthetic full-axis mean-contract tests before biological scoring. The executed sampler combines multinomial draws,16 depth-stratum gene-total repairs and bounded integer swaps. All output cell libraries and reference gene totals are exact; mean per-cell expression drift is measured, not assumed zero. Repairs can alter shape, so diagnostics use the final integer counts. A successful mean contract verifies isolation of the mean predictor; it does not establish biological realism.

## Independent audit

- Four focused neural/count tests passed before the plan was frozen.
- {audit['source_count_contexts_checked']} source contexts and{audit['original_full_axis_rows_checked']} complete original cell rows were checked on the full gene axis.
- {audit['training_only_group_means_checked']} fitting-only group means were independently recomputed; validation cells never supplied these means.
- Training/validation identities, state-shuffle batch/partition membership, family exclusions, checkpoint hashes and finite epoch traces passed.
- All executed null-cache hashes, control labels and disjoint actual/reference/inference identities passed.
- {audit['full_axis_count_replays']} full-axis baseline/conditional count replays and their scorer call counts matched exactly.
- Frozen calibration checks were independently recomputed; the automatic stop follows the saved failed gate.

No treatment comparison, later-family training, reserve score or host score is claimed. The pipeline stops at the first failed required calibration as specified before fitting. Existing unexecuted code paths are documented as prepared, not validated. The remaining study/family limits, observational state proxies and source-to-query distribution shift remain relevant.

## Where the remaining errors occur

A separate exact replay of the conditional null locates its352 called genes:

|Coordinate category|All genes|Eligible genes|Null DE calls|
|---|---:|---:|---:|
|Trained; >=50 source detections|15215|10209|142|
|Trained; <50 source detections|2257|89|4|
|Panel genes with empirical fallback|300|292|109|
|Other empirical fallback|761|191|97|

Thus206/352 calls (58.5%) lie in1061/18533 fallback coordinates (5.7% of the output axis). The trained coordinates still contribute146 calls, so changing the fallback alone is not demonstrated to pass calibration. These categories involve different genes; this is descriptive localization, not an isolated causal ablation proving which operation caused each error. The null-gene diagnostic matches the saved conditional count digest exactly. Per-gene q-values, LFCs, detection, variance and source support are in H1_null_called_genes.csv; category counts are in H1_null_error_categories.csv.

This creates two concrete next hypotheses: avoid adding another sampling layer to already observed fallback counts, and adapt learned control distributions using the query context's own available NTCs. Panel genes need exclusion from leakage-prone state features; masking all300 from learned output variation is a conservative choice tested here, not a competition requirement. Any revised mask, sampler or adaptation procedure needs a new plan and validation, rather than changing this failed V26 after seeing outcomes.

## Reporting recovery and latent check

A Path/string mismatch interrupted saving the first completed assessment. The operational assess_v26_recovery.py adapter coerces hash paths, executes the unchanged hash-pinned assessor and preserves table hashes. Original code, fitted arrays, populations and thresholds remain unchanged. run_v26_with_recovery.py is the from-scratch entry point; after recovery the frozen pipeline resumed verified caches and recorded its planned stop. The independent auditor separately recomputes every failed check.

Descriptive source-control diagnostics use1276 validation cells. Conditional aggregated posterior variance ranges0.919..1.062 around the normal generation prior's1, with mean-vector norm0.102. This does not indicate a large diagonal source-prior variance collapse; it does not verify full joint distributions or H1 transfer. See H1_latent_diagnostic.json and root diagnose_v26_latents.py. No diagnostic changes the release rule.

## Implication

This experiment directly tests a complete trainable count-distribution mechanism instead of adding bounded residuals to the old generator. Better source likelihood cannot substitute for calibrated unseen-context distributions. This configuration has not earned a treatment or leaderboard release. Any adaptation to the query controls, change to the likelihood/latent prior or mean-repair rule must be a separately frozen experiment, preserving this result. The present evidence does not identify a hidden competition failure or promise improvement from another architecture. A justified next experiment would calibrate the generator to each query context's available NTC distributions, then learn source-supported treatment-induced distribution changes. This needs a new frozen protocol and matched controls.

## Reproduce

Use the unchanged scientific Python3.12 runtime and versions in analysis/runtime.json. Existing frozen plans refuse overwrite. Source/weight/count arrays are regenerable and ignored by Git. Root scripts and small evidence are versioned. Run:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
set "VCC_SCI_PY=C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe"
"%VCC_SCI_PY%" audit_v26_experiment.py
"%VCC_SCI_PY%" diagnose_v26_latents.py --family H1
"%VCC_SCI_PY%" diagnose_v26_null_genes.py
"%VCC_SCI_PY%" document_v26_results.py
```

From scratch, V26_STRATEGY.md lists preparation and gated training commands. Replaying the saved pipeline resumes identical fitted weights/caches and preserves its stop. It does not revise parameters or reopen treatment outcomes.
''';(ROOT/'V26_RESULTS.md').write_text(report,encoding='utf-8');(ROOT/'V26_SUBMISSION.md').write_text(f'''# V26 artifact decision

No V26 submission is released. The source-trained conditional count decoder fails frozen real-control calibration in{key}, so the pipeline stops before treatment scoring or reserve evaluation. See V26_RESULTS.md for actual measurements and V26_STRATEGY.md for methods.

Retain the already published V24 entry `{champion['entry_id']}`, overall **{champion['overall']}**. Its package remains `experiments/v24_cd4_response/prediction.vcc`. A repeated upload would duplicate the current reference. No V26 upload, duplicate champion copy or Git push is performed.

Source data, fitted neural arrays and count caches remain local/ignored. The root code, README, frozen plan, model traces, diagnostic tables and independent audit are committed for reproduction. Previous numerical files/packages are unchanged.
''',encoding='utf-8');section=f'''<!-- v26-count-decoder:begin -->

## V26 neural count-distribution test (2026-10-15)

Implemented a conditional logistic-normal multinomial decoder with fixed V24 response means and exact cell libraries/gene totals. {len(models)} networks were fitted with{training['rows']:,} sampled training/validation source cells under complete{key} exclusion. The executed real-control null calibration fails its frozen gate; treatment outcomes and the remaining donor reserve stay unopened. Mean constraints and numerical/input/count replay audits are recorded separately from statistical usefulness. No new package is released; V24 remains best recorded at{champion['overall']}.

Root implementation and reproducibility: [V26 strategy](V26_STRATEGY.md). Actual diagnostics and limitations: [V26 results](V26_RESULTS.md). [Artifact decision](V26_SUBMISSION.md). Prior numerical models/packages are preserved; large source arrays and weights stay outside Git.

<!-- v26-count-decoder:end -->
''';readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8');a='<!-- v26-count-decoder:begin -->';z='<!-- v26-count-decoder:end -->'
 if a in text:text=text[:text.index(a)]+section.strip()+text[text.index(z)+len(z):]
 else:
  lines=text.splitlines(keepends=True);text=lines[0]+'\n'+section+'\n'+''.join(lines[1:])
 readme.write_text(text,encoding='utf-8');print('V26 measured calibration rejection documented; V24 retained',flush=True)
if __name__=='__main__':main()
