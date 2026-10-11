"""Render measured V15 research and upload documentation from verified evidence."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v15_weighted_transport/analysis';OUT=A.parent

def table(frame):
    def val(x):
        if isinstance(x,(float,np.floating)):return 'NA' if not np.isfinite(x) else f'{x:.6g}'
        return str(x)
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+['| '+' | '.join(val(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])

def read(name):return json.loads((A/name).read_text())
def main():
    assert read('development_evaluation.json')['complete'] and read('null_evaluation.json')['complete'];selection=read('selection.json');audit=read('independent_experiment_audit.json');assert audit['passed'] and audit['original_baseline_same_v15_input_count_and_DE_replays']==18 and len(audit['count_profile_caches'])==34
    ready=read('submission_ready.json') if (A/'submission_ready.json').exists() else None
    if ready:assert read('independent_validation.json')['passed'] and ready['sha256']==build.sha256(OUT/'prediction.vcc')
    confirmation=read('confirmation_result.json') if (A/'confirmation_result.json').exists() else None
    matched=read('matched_control_result.json') if (A/'matched_control_result.json').exists() else None
    summary=pd.read_csv(A/'development_summary.csv');panels=pd.read_csv(A/'development_panels.csv')
    source=pd.read_csv(A/'response_summary.csv');div=pd.read_csv(A/'expert_diversity.csv')
    diversity=div.groupby(['family','expert'])[['error_correlation','corrects_v13_sign','breaks_v13_sign']].mean().reset_index()
    models=[]
    for file in sorted((A/'models').glob('*.json')):
        m=json.loads(file.read_text())['model'];models.append({'held_family':file.stem,'V13_weight':m['weights'][0],'Student_weight':m['weights'][1],'specialist_weight':m['weights'][2]})
    comparisons=pd.DataFrame([{'candidate':r['candidate'],'passes_available_guards':r['pass'],'complete_primary':r['complete_primary_metrics'],**r['changes']} for r in selection['comparisons']])
    counts={'null':read('null_evaluation.json')['population_cases'],'development':read('development_evaluation.json')['population_cases']}
    if (A/'confirmation_evaluation.json').exists():counts['reserved_verification']=read('confirmation_evaluation.json')['population_cases']
    if matched and matched.get('applicable'):counts['supplementary_matched_control']=matched['generated_control_populations']
    status='A new package passed official prep and the independent full-payload audit. It has not been uploaded and has no leaderboard score.' if ready else 'No new package was released; V13 remains the upload reference. The completed comparisons did not justify replacing it under the frozen release policy.'
    response_comparison=next(r for r in selection['comparisons'] if r['candidate']=='response')
    changes=response_comparison['changes']
    decision=(f"The compact response blend changed fidelity by {changes['fidelity']:+.6f}, reach by {changes['reach']:+.6f}, Jaccard by {changes['jaccard']:+.6f}, and eligible NMAE by {changes['nmae']:+.6f}. Its paired expression error worsened by {100*changes['relative_expression_worsening']:.2f}%. These are family-balanced public changes; lower NMAE is better.")
    confirmation_text=table(pd.read_csv(A/'confirmation_summary.csv')) if confirmation else 'Reserved verification was not opened after the earlier rejection.'
    verified='passed' if confirmation and confirmation['new_candidate_passes'] else 'not passed / not reached'
    text=f'''# V15 results: weighted response and population transport

{status}

V13 remains the best *published* project score, 0.137883134. Public raw metrics below do not predict an official leaderboard score or a rank gain. Numerical V13/V14 dependencies and prior packages were preserved.

{decision}

The tested ensemble did not supply enough independent error correction. The distribution variants increased apparent DE calls in control-only comparisons, and their response-count improvements did not survive the full guard policy. No reserved verification or production fit is opened after a failed development gate.

## Implemented experiment

Three response experts were tested: exact V13, Student-t pooling with neutral unknown-knockdown handling, and a training-source-reliability/similarity specialist. A six-point convex grid retains >=50% V13, with weights calibrated through inner-family holdouts inside each outer-held-family evaluation. This is weighted continuous-effect prediction, not a majority vote on counts. Missing measurements abstain; CD4 conditions remain one family. Per-gene direction agreement is diagnostic evidence, not independent model votes.

Source-treated/control marginal changes were extracted from verified H1/HCT116/HEK293T counts. Five relative quantiles in two depth bins are transferred onto ranked empirical query-control cells, then blended at 10% or 25% on measured supported genes and integer-fitted to the same response moments. This tests one gene-specific response-distribution construction. The linear blend also attenuates the original Poisson term by (1 - blend): its coefficient falls from 0.65 to 0.585 or 0.4875 on supported genes before fitting. It therefore changes existing sampling noise as well as adding empirical variation. The zero-response null uses zero quantile deltas and still fails, so these results do not isolate a failure of cross-source quantile transfer itself. A future noise-preserving construction would need a new frozen experiment, not a retrospective change to V15. Extreme quantiles, pooled source batches, assay differences and finite treatment counts remain limitations. See [V15_STRATEGY.md](V15_STRATEGY.md) for the exact equations, thresholds, identities and commands.

Population-case counts: **{counts}**. These include repeated seeds, multiple predictions for shared truths, and null/control diagnostics; they are not that many independent biological experiments. Development uses 14 targets with 400 genuine cells per case; reserved algorithm verification uses nine other current target names. Older public biology was already inspected, so this is not blind biological confirmation.

## Ensemble screen and complementarity

Frozen response for count testing: **{read('response_selection.json')['response']}**.

{table(source)}

Nested mixture weights:

{table(pd.DataFrame(models))}

Error correlations and direction corrections/breaks are averaged per source family/expert. The direction counts describe reference-DE genes per target/context contrast, not all-gene accuracy or official fidelity. High residual correlation limits the diversification benefit of an ensemble.

{table(diversity)}

All response-screen metrics are public summary proxies. Whole-family exclusion prevents using held-family intervention labels in a predictor; previously inspected source biology still limits inferential novelty. Source-specialist weights of zero mean that expert was rejected, not secretly included as an independent vote.

## Disjoint NTC diagnostics

{table(pd.read_csv(A/'null_summary.csv'))}

Predicted calls compare generated controls with independent right-hand controls. The actual NTC diagnostic cells are disjoint from templates but are included in the left query-mean fitting pool; their expression numerator is optimistic and is not the null selection threshold. A distribution finalist may not increase generated-versus-independent-control DE calls by >10% relative to compact emission. NTC diagnostics use 200 cells per arm, not a 400-cell challenge treatment test.

## Actual-count development

Higher is better for PDS/fidelity/reach/Jaccard; lower for NMAE/expression ratio. Family-balanced NMAE averages only matched eligible panels. Large expression ratios can be dominated by tiny control-distance denominators; the predeclared guard instead compares paired absolute numerators, equivalent to relative ratio worsening when the shared denominator resolves. Missing primary metrics are disclosed rather than filled with zero.

{table(summary)}

Available-metric guard comparisons:

{table(comparisons)}

Frozen finalist: **{selection['candidate']}**. Equal metric-improvement ranks select among passing candidates; these ranks are not private competition anchor scaling. Full check details and missing-NMAE panels are in `selection.json`. The 2% expression tolerance is an exploratory prospective bound, not statistical noninferiority and not a retroactive change to V14's outcome.

## Matched negative control and reserved verification

Matched distribution control: **{matched if matched else 'not reached'}**.

The original shuffled25 development arm matches response/quantile25 exactly, but differs in blend from response/quantile10 and in means from V13/quantile25. Before generated-count outcomes, an extra freeze required an exactly mean/blend-matched shuffled-target check for any distribution finalist. It can reject the finalist without retuning or selecting another from the same outcome.

Reserved algorithm verification: **{verified}**.

{confirmation_text}

A pass on available guards is not complete six-metric verification when NMAE panels are missing. Current target-name splits and paired seeds do not create new biological lineages. No candidate is reselected from a failed reserved result.

## Integrity and reproduction

Nine scientific tests passed (five new V15 tests plus four imported V14 invariants). The independent experiment audit verifies original unique IDs and disjoint splits, manually reconstructs **{sum(r['quantile_coordinates'] for r in audit['quantile_replays'])}** source quantile coordinates, checks nested grid choices/cache hashes, and exactly replays **{audit['original_baseline_same_v15_input_count_and_DE_replays']}** original-baseline count/profile/DE cases on V15's exact inputs. The initial historical-hash assertion failed: V14 sampled unsorted control pools before saving sorted pools, whereas V15 sampled the saved pools. Actual treatment identities and whole control halves match, but template identities differ. The failed auditor is preserved as `audit_v15_experiment_initial.py`; `audit_correction.json`, `historical_baseline_differences.csv` and `same_input_baseline_replay.json` record the correction. Frozen methods, V15 identities, outcomes and candidate selection were not changed. This separates numerical integrity from a performance claim. If a package is present, its independent release audit also checks every payload count/index/pointer/axis, every sampled original control identity/guide/library, all production response/distribution arrays and twelve complete-group regenerations.

Root scripts: `v15_ensemble_methods.py`, `v15_population_methods.py`, `prepare_v15_experiment.py`, `run_v15_experiment.py`, `assess_v15_experiment.py`, `check_v15_matched_control.py`, `audit_v15_experiment.py`, `replay_v15_baseline.py`, `audit_v15_experiment_initial.py` (preserved failed audit), `fit_v15_submission.py`, `build_v15_submission.py`, `package_v15_submission.py`, `validate_v15_submission.py`, `release_v15_submission.py`, `test_v15_ensemble.py` and this report generator. Exact commands are in [V15_STRATEGY.md](V15_STRATEGY.md). The release scheduler only sequences frozen stages and stops on rejection; it performs no upload.

Plans, source/model sidecars, scores, selections and audits are under `experiments/v15_weighted_transport/analysis/`. Large arrays, generated counts, packages and full sampling provenance are regenerated locally and excluded from Git. `report_manifest.json` pins the evidence and generated report bytes.
'''
    result=ROOT/'V15_RESULTS.md';result.write_text(text,encoding='utf-8',newline='\n')
    outputs={'V15_RESULTS.md':build.sha256(result)}
    if ready:
        model=read('reference_model.json');validation=read('independent_validation.json');name='Robust response ensemble v15' if model['emitter']=='compact' else 'Robust ensemble and source-linked transport v15'
        submit=f'''# V15 upload-ready submission

Ready at `experiments\\v15_weighted_transport\\prediction.vcc`. Official VCC prep and the independent complete-payload audit passed. **Exploratory prediction; no leaderboard result yet.** V13 remains the best published reference. See [V15_RESULTS.md](V15_RESULTS.md) for public gains, failed/incomplete checks and limitations.

- Candidate: `{model['candidate']}`; response `{model['response']}`; emitter `{model['emitter']}`.
- Size: {ready['bytes']:,} bytes ({ready['bytes']/2**30:.2f} GiB).
- SHA-256: `{ready['sha256']}`.
- Shape: 360,000 cells x 18,533 genes; 900 context/target groups, exactly 400 cells per group.
- Every payload count/index/pointer/axis compared; all sampled original control identities and library totals verified; {validation['exact_replayed_groups']} complete groups regenerated exactly.
- Complete six-metric public verification: {model['confirmed_improvement']}. Available-public-guard verification: {model['available_public_guards_passed']}. Missing/ineligible NMAE panels and reused biology prevent stronger conclusions.

Run these commands in **Command Prompt**, which matches your earlier terminal:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v15_weighted_transport\\prediction.vcc -m "{name}" --wait
```

If the terminal is interrupted after upload, check the printed entry ID rather than submitting the same file again:

```cmd
.\\.competition_env\\Scripts\\python.exe -m vcc status ENTRY_ID --wait
```

The file has not been uploaded by this workflow. Do not interpret format/integrity checks as a guaranteed score or rank improvement.
'''
        path=ROOT/'V15_SUBMISSION.md';path.write_text(submit,encoding='utf-8',newline='\n');outputs[path.name]=build.sha256(path)
    else:
        baseline=read('baseline_package_check.json');assert baseline['passed']
        path=ROOT/'V15_SUBMISSION.md'
        retained=f'''# V15 release decision: retain V13

No new V15 package was created. All four deployable development candidates failed the frozen release criteria; see [V15_RESULTS.md](V15_RESULTS.md) for measured comparisons and [V15_STRATEGY.md](V15_STRATEGY.md) for reproduction. Reserved verification was not opened and no failed candidate was reselected.

The retained file is `experiments\\v13_response_programs\\prediction.vcc`, the existing V13 submission already published as `AvAsYkuxFSjbxYDCd0wJ` with overall 0.13788313368012328. Its SHA and official format were checked again after the V15 tests. **There is no reason to re-upload this identical package.**

- Size: {baseline['bytes']:,} bytes ({baseline['bytes']/2**30:.2f} GiB).
- SHA-256: `{baseline['sha256']}`.
- Shape: 360,000 cells x 18,533 genes; 900 groups of 400 cells.
- Artifact check: `experiments/v15_weighted_transport/analysis/baseline_package_check.json`.

If this file is needed for a different permitted evaluation, its Command Prompt command is:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v13_response_programs\\prediction.vcc -m "Assay-aware knockdown transport v13" --wait
```

No upload was performed by V15. The retained V13 model was exploratory when built; its published score is the observed evidence, not a guarantee on the final test set. A format check does not validate biological accuracy.
'''
        path.write_text(retained,encoding='utf-8',newline='\n');outputs[path.name]=build.sha256(path)
    marker='<!-- v15-research:begin -->';end='<!-- v15-research:end -->'
    outcome='A new exploratory package passed official prep and independent complete-payload validation; see [V15_SUBMISSION.md](V15_SUBMISSION.md) for the Command Prompt upload command.' if ready else 'No new package was released; V13 remains the upload reference. See [V15_SUBMISSION.md](V15_SUBMISSION.md) for the retained artifact and why no repeat upload is needed.'
    block=f'''{marker}

## V15 response ensemble and population transport

Implemented nested held-family calibration of V13, robust Student-t pooling and a source specialist, plus depth-stratified source-treated marginal changes and matched shuffled controls. Frozen development includes 168 generated populations from 14 targets/two seeds; separate NTC diagnostics and reserved algorithm verification guard release. Selected development candidate: `{selection['candidate']}`. {outcome} These are public exploratory comparisons on partly reused biology, with missing-metric and tiny-denominator limitations; no new leaderboard score is claimed.

Methods/reproduction: [V15_STRATEGY.md](V15_STRATEGY.md). Measured comparisons, expert diversity and verification: [V15_RESULTS.md](V15_RESULTS.md). Root code and small hashed evidence are versioned; prior numerical models/packages remain unchanged.

{end}'''
    readme=ROOT/'README.md';s=readme.read_text(encoding='utf-8')
    if marker in s:
        i=s.index(marker);j=s.index(end,i)+len(end);s=s[:i]+block+s[j:]
    else:s=s.replace('<!-- v14-research:end -->','<!-- v14-research:end -->\n\n'+block,1)
    readme.write_text(s,encoding='utf-8',newline='\n')
    names=['audit_correction.json','same_input_baseline_replay.json','historical_baseline_differences.csv','baseline_package_check.json','plan.json','matched_control_plan.json','unit_tests.json','response_selection.json','response_summary.csv','expert_diversity.csv','null_evaluation.json','null_summary.csv','development_evaluation.json','development_scores.csv','development_panels.csv','development_summary.csv','selection.json','independent_experiment_audit.json','matched_control_result.json','matched_control_scores.csv','confirmation_result.json','confirmation_summary.csv','reference_model.json','submission_ready.json','independent_validation.json']
    build.write_json(A/'report_manifest.json',{'inputs':{str((A/n).relative_to(ROOT)):build.sha256(A/n) for n in names if (A/n).exists()},'outputs':outputs,'generator_sha256':build.sha256(Path(__file__)),'uploaded':False,'scope':'Measured public-source research and, if present, independently verified exploratory package; no private score or rank guarantee.'})
    print('V15 results/README/documentation generated; package present:',bool(ready),flush=True)
if __name__=='__main__':main()
