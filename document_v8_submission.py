"""Write V8 outcome and reproduction documentation after verified packaging."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import build_submission as build
from v8_methods import ROOT,A,OUT


def table(frame):
    def value(x):
        if isinstance(x,(float,np.floating)):return f'{x:.6f}' if np.isfinite(x) else 'undefined'
        return str(x).replace('|','/')
    return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |',
                      *['| '+' | '.join(value(x) for x in row)+' |' for row in frame.itertuples(index=False,name=None)]])


def main():
    ready=json.loads((A/'submission_ready.json').read_text());audit=json.loads((A/'independent_validation.json').read_text())
    assert audit['passed'] and audit['package_sha256']==ready['sha256']
    model=json.loads((A/'reference_model.json').read_text());safe=json.loads((A/'safe_count_selection.json').read_text());confirmation=json.loads((A/'safe_confirmation_result.json').read_text())
    plan=json.loads((A/'validation_plan.json').read_text());neutral=pd.read_csv(A/'decoder_neutral_summary.csv')
    cal=pd.read_csv(A/'count_calibration_summary.csv');conf=pd.read_csv(A/'count_confirmation_summary.csv')
    fields=['candidate','pds','fidelity','reach','jaccard','nmae','common_expression_ratio','utility']
    shape=ready['prediction_validation'];effect=A/'reference_model_arrays.npz'
    comparison=None
    v7=ROOT/'experiments/v7_multisource_adaptive/analysis/reference_model_arrays.npz'
    if v7.exists():
        with np.load(effect,allow_pickle=False) as z,np.load(v7,allow_pickle=False) as prior:
            difference=float(np.max(abs(z['effects']-prior['effects'])))
            comparison={'previous_v7_response_sha256':build.sha256(v7),'current_effects_sha256':build.sha256(effect),
                        'maximum_absolute_log2_difference':difference,'exactly_same_response_effects':bool(difference==0)}
        build.write_json(A/'response_comparison.json',comparison)
    passed=confirmation['confirmed_improvement']
    contexts=[]
    for c in plan['contexts']:
        cells=c['cells'];contexts.append({'context':c['context'],'left_controls':len(c['left']),'right_controls':len(c['right']),
            'calibration_targets':len(set(v['target'] for v in cells if v['split']=='calibration')),
            'confirmation_targets':len(set(v['target'] for v in cells if v['split']=='confirmation')),
            'observed_cells_min':min(len(v['actual']) for v in cells),'observed_cells_max':max(len(v['actual']) for v in cells),
            'old_v7_count_targets_reused':len(c['old_v7_targets_reused'])})
    pd.DataFrame(contexts).to_csv(A/'validation_coverage.csv',index=False)
    availability=[]
    for stage in ['calibration','confirmation']:
        panels=pd.read_csv(A/f'count_{stage}_panels.csv')
        for (family,candidate),rows in panels.groupby(['family','candidate']):
            availability.append({'stage':stage,'family':family,'candidate':candidate,'panels':len(rows),
                'finite_utility_panels':int(rows.utility.notna().sum()),
                'resolved_expression_panels':int(rows.common_expression_ratio.notna().sum()),
                'finite_reach_panels':int(rows.reach.notna().sum())})
    availability=pd.DataFrame(availability);availability.to_csv(A/'metric_availability.csv',index=False)
    text=f'''# V8 submission: shared-control validation and calibrated count profiles

Upload-ready file: `experiments/v8_corrected_transfer/prediction.vcc` ({ready['bytes']:,} bytes). It contains **360,000 cells x 18,533 genes**, all **300 targets x 400 cells x A/B/C**. Official VCC prep and independent full-payload validation passed. This package has not been uploaded and has no leaderboard score.

SHA-256: `{ready['sha256']}`.

## Final decision

Candidate **{model['candidate']}**, response **{model['response']}**, count generator **{model['generator']}**. The response profile screen favored `v4_a65`; the generated-cell safety selection chose the candidate above. Confirmation improvement under the declared public proxy: **{passed}**, utility difference **{confirmation['utility_difference']:+.6f}**. Confirmation was assessed once without reselection. Public proxy utility is not an official anchored leaderboard score.

The original count-only calibration choice `{safe['original_count_winner']}` was not accepted automatically. An independent neutral-control test rejected `dual_half`: it produced 968-2,403 significant genes without any applied perturbation. Both `legacy` and `dual` produced only 0-1 false calls across the six control tests. V8 uses the safe frozen calibration winner.

{('The response-effect arrays are exactly equal to the V4/V7 fallback arrays (maximum difference 0). This is an isolated count-decoder change, not a new biological response model. Unlike the V7 resampling comparison, it changes the expected mean/pooled profiles and repairs both integer count moments.' if comparison and comparison['exactly_same_response_effects'] else 'The response-effect arrays differ from the V7 fallback; the response and decoder choices are recorded separately.')}

## What was fixed

- One common reference control profile and cell group per evaluation stratum/seed replaces target-specific control origins. Constant predictions now score exactly 0.5.
- Control-only batch PCA/KMeans defines technical strata before response modeling. Prediction templates use a target-independent 400-cell draw, never actual treated batch proportions.
- Source lineages are excluded from their own response training. Calibration and confirmation target genes are globally disjoint. None of these selected targets were used in V7 generated-cell tests; they remain public development data already available to source-summary work.
- Expression calculations retain signed corrections, scalar whole-gene jackknifes, and explicit unresolved-denominator reporting. DE uses the same common control, its CPM>5 gene universe, and per-target BH correction.
- Response pooling, shrinkage and amplitude were screened separately from the count decoder. Full moment matching anchors mean cell composition and pooled expression to the entire available control population, preserves cell depths, and verifies realized integer moments.
- Neutral/control and independent-replicate diagnostics are recorded. A reduced-variance decoder was rejected for creating extensive false DE rather than adopted for its apparent metric gain.

## Public calibration

{table(cal[fields])}

## Public confirmation

{table(conf[fields])}

Every metric is shown, including tradeoffs. Fidelity is the raw yield metric, not sign precision. Many HCT116/HEK293T development targets have zero or tiny reference DE sets, which makes unweighted fidelity comparisons and candidate-specific undefined values unstable. Before confirmation, the final safety selector retained expression/PDS guards, required every neutral false-DE fraction <=0.5%, and replaced its initial mean-fidelity guard with a common reference-supported budget-weighted fidelity guard (reference calls >=10; within 0.02 of baseline). The original calibration decision and all raw metric means remain preserved. This is a calibration refinement, not a claim that the budget-weighted guard is the official leaderboard metric.

## Metric availability

{table(availability)}

Undefined values remain undefined. Aggregation averages defined panels/strata/families separately for each metric; if a family's utility is entirely undefined, its missing value is skipped. A nominal equal-family proxy can therefore include fewer than three lineages. The raw CSVs and availability table expose those omissions. A positive corrected expression denominator and a nonempty real DE budget are required for the corresponding ratio/reach metrics.

## Neutral decoder checks

{table(neutral)}

## Validation coverage and limits

{table(pd.DataFrame(contexts))}

The source tests contain 60-200 actual cells per target, with 400 generated and reference cells. Technical strata reduce batch variation but cannot establish matching platform, timing, knockdown efficiency or response biology for hidden A/B/C treatments. Replicate diagnostics use independent treated halves and independent control origins; their smaller group sizes are not official anchors. Public scores, including Jaccard gains on empty-response targets, must not be presented as hidden performance estimates. {('The frozen safe candidate improves the declared confirmation proxy; a leaderboard gain remains unknown.' if passed else 'The frozen safe candidate did not satisfy the confirmation improvement rule. This upload package is an experimental calibration candidate, not a demonstrated improvement over V4.')}

## Reproduce

From Windows Command Prompt, with existing pinned source/control caches:

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe run_v8_pipeline.py --resume --workers 4
```

For a fresh reproduction, omit `--resume` in a checkout without V8 outputs. Stages are: unit tests, control-only plan, response calibration, generated-cell calibration, neutral decoder checks, safety selection, confirmation, final confirmation assessment, production fit, model audit, preflight, full count generation, official prep, independent validation, documentation. Original V4-V7 numerical modules and payloads are preserved. A fresh numerical rerun requires a workspace without the V8 experiment directory; archive the recorded V8 directory before starting rather than overwriting its frozen records.

Key root scripts: `prepare_v8_validation.py`, `v8_methods.py`, `evaluate_v8_submission.py`, `benchmark_v8_neutral.py`, `select_v8_safe_generator.py`, `confirm_v8_safe_selection.py`, `fit_v8_submission.py`, `build_v8_submission.py`, `validate_v8_submission.py`, `test_v8_submission.py`, `run_v8_pipeline.py`, `document_v8_submission.py`.

## Submit

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v8_corrected_transfer\\prediction.vcc -m "Shared-control calibrated count profiles v8" --wait
```

If interrupted after upload/scoring starts, use `vcc status <entry-id> --wait` instead of creating another entry.

## Integrity records

`submission_ready.json` holds the package checksum and official prep metadata. `independent_validation.json` confirms every packaged count/index/pointer and axis matches H5AD, scans all 360,000 cell libraries, checks all control IDs and guide sampling, and exactly replays 12 groups (88,958,400 gene-count entries). `reference_model.json` freezes model arrays, code and inputs; selections and evaluation files preserve calibration/confirmation chronology and sampling. `analysis/v8_official_metric_review.json` pins the official source examined for metric parity.
'''
    (ROOT/'V8_SUBMISSION.md').write_text(text,encoding='utf-8')
    readme=ROOT/'README.md';content=readme.read_text(encoding='utf-8-sig')
    start='<!-- v8-status:begin -->';end='<!-- v8-status:end -->'
    block=f'{start}\n\nV8 is retained as a verified diagnostic pilot at `experiments/v8_corrected_transfer/prediction.vcc`, with official prep and independent payload verification passed. Candidate `{model["candidate"]}` uses `{model["response"]}` and `{model["generator"]}`; public confirmation proxy improvement accepted: {passed}. It has not been submitted. See [V8_SUBMISSION.md](V8_SUBMISSION.md) for the corrected validation, neutral controls, tradeoffs and upload command.\n\n{end}'
    if start in content:content=content[:content.index(start)]+block+content[content.index(end)+len(end):]
    else:content=content.replace('## Completed work',block+'\n\n## Completed work',1)
    readme.write_text(content,encoding='utf-8')
    print('V8 submission documentation and README updated',flush=True)

if __name__=='__main__':main()
