"""Document verified V9 artifacts and all calibration/confirmation tradeoffs."""
import json
import pandas as pd
import build_submission as build
from v9_methods import ROOT,A,OUT
from document_v8_submission import table


def main():
    ready=json.loads((A/'submission_ready.json').read_text());audit=json.loads((A/'independent_validation.json').read_text())
    assert audit['passed'] and audit['package_sha256']==ready['sha256']
    model=json.loads((A/'reference_model.json').read_text());selection=json.loads((A/'count_selection.json').read_text())
    result=json.loads((A/'confirmation_result.json').read_text());plan=json.loads((A/'validation_plan.json').read_text())
    cal=pd.read_csv(A/'count_calibration_summary.csv');conf=pd.read_csv(A/'count_confirmation_summary.csv')
    stable_fidelity=pd.DataFrame([{'stage':stage,'candidate':candidate,'supported_budget_fidelity':value}
        for stage,record in [('calibration',selection),('confirmation',result)]
        for candidate,value in record['supported_budget_fidelity'].items()])
    stable_fidelity.to_csv(A/'supported_fidelity_summary.csv',index=False)
    neutral=pd.read_csv(A/'decoder_neutral_summary.csv');fields=['candidate','pds','fidelity','reach','jaccard','nmae','common_expression_ratio','utility']
    coverage=[]
    for c in plan['contexts']:
        coverage.append({'context':c['context'],'left_controls':len(c['left']),'right_controls':len(c['right']),
            'calibration_targets':len({v['target'] for v in c['cells'] if v['split']=='calibration'}),
            'confirmation_targets':len({v['target'] for v in c['cells'] if v['split']=='confirmation'}),
            'actual_cells_min':min(len(v['actual']) for v in c['cells']),'actual_cells_max':max(len(v['actual']) for v in c['cells'])})
    pd.DataFrame(coverage).to_csv(A/'validation_coverage.csv',index=False)
    availability=[]
    for stage in ['calibration','confirmation']:
        panels=pd.read_csv(A/f'count_{stage}_panels.csv')
        for (family,candidate),rows in panels.groupby(['family','candidate']):
            availability.append({'stage':stage,'family':family,'candidate':candidate,'panels':len(rows),
                'finite_utility_panels':int(rows.utility.notna().sum()),
                'resolved_expression_panels':int(rows.common_expression_ratio.notna().sum()),
                'finite_reach_panels':int(rows.reach.notna().sum())})
    availability=pd.DataFrame(availability)
    availability.to_csv(A/'metric_availability.csv',index=False)
    command=r'.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v9_activation_preserving\prediction.vcc -m "Activation-preserving profile calibration v9" --wait'
    text=fr"""# V9: activation-preserving profile calibration

Upload-ready: `experiments/v9_activation_preserving/prediction.vcc` ({ready['bytes']:,} bytes).
**360,000 cells x 18,533 genes; 300 targets x 400 cells x three contexts A/B/C.**
Official VCC prep and independent full-payload verification passed. No upload or leaderboard score yet.
SHA-256: `{ready['sha256']}`.

## Decision and evidence

Frozen generator: **{model['generator']}**; biological effects: **exact unchanged V4/V7 response arrays**.
Calibration chose the best six-metric public proxy among candidates passing the declared neutral, fidelity, PDS and expression guards. Calibration froze before confirmation. Confirmation improvement: **{result['confirmed_improvement']}**, guard passed **{result['guard_passed']}**, proxy utility difference **{result['utility_difference']:+.6f}**. No confirmation-based reselection occurred. These are unanchored public metrics, not a predicted leaderboard overall.

## Why this experiment

The prior validation audit found target-specific control origins inflated constant-profile PDS above its correct 0.5 null. Expression corrections also needed whole-gene scalar jackknifes and signed denominators. V8 fixed those calculations, used control-defined technical strata, target-independent templates, source-lineage exclusion and globally disjoint target calibration/confirmation. Its safe full moment decoder improved some expression metrics but lost supported-target fidelity and did not confirm a proxy improvement. A half-variance decoder produced 968-2,403 false DE calls on unperturbed controls and was rejected. See `V8_SUBMISSION.md` and `TRANSFER_VALIDATION_AUDIT.md`.

V9 retains these validation repairs and tests a different count initialization. First simulate V4's binomial downregulation and zero-aware Poisson activation. Then calibrate mean cell composition and pooled expression by 25%, 50% or 100% toward the full control-population response. Calibration starts from the simulated perturbed cells, preserving activation of genes that were zero in the original template. It preserves each simulated cell's library exactly, rather than imposing untreated depths. Integer rounding and bulk repair are upstream reviewed implementations, with verified realized moment errors. Legacy generation remains an explicit baseline.

The correction strength is selected on new count-validation target genes, excluding every V8 tested target globally and V7 count-tested targets in each held source. All targets were already public source-summary data; these are new generated-cell tests, not pristine blinded biology. HCT116/HEK293T eligibility requires >=150 actual cells in the control-defined stratum; H1 requires >=60. Selection uses metadata and own control expression >5 CPM, never treated expression values. Fixed 400-cell template and reference groups are shared across all targets per context/seed, from disjoint untreated populations. No treated batch proportions enter model inputs.

H1 evaluation excludes direct H1 components and its learned prior; HCT/HEK are absent from V4's response. K562/CD4 supply the transferred fallback. A/B/C still have no measured treated training labels. Consequently the central cross-lineage causal transfer uncertainty remains unresolved; control gene reconstruction neural networks cannot supply those missing intervention responses.

## Calibration: all metrics

{table(cal[fields])}

## Confirmation: all metrics

{table(conf[fields])}

Raw fidelity follows the official direction-yield budget, not sign precision. Undefined zero-budget values remain undefined. The safety guard uses the same reference-call>=10 targets for each candidate and count-weighted fidelity within 0.02 of legacy; expression <=1.02 legacy, PDS within 0.05, and every neutral-control false-DE fraction <=0.005. This stability guard is a development criterion, not a replacement official metric. Utility is the V8 equal-family average of centered PDS, fidelity, reach, Jaccard, clipped 1-NMAE and clipped 1-expression. Raw expression denominators can be weak, and absent hidden anchor sets prevent translating these values into leaderboard scores. Reference-supported fidelity and guards are recorded separately in the JSON decisions.

## Reference-supported fidelity guard

{table(stable_fidelity)}

Passing a tolerance guard is not strong absolute causal accuracy. Confirmation supported-target fidelity is very low for both legacy and calibrated transfer; this remains a major response-model weakness. The improved decoder does not supply missing matched A/B/C intervention data or establish that the V4 response directions are correct. The direct target retains V4's fixed 90% knockdown assumption; it is excluded from scored gene sets but still affects library normalization. This experiment does not learn hidden dose or timing. Raw mean fidelity can be dominated by tiny reference budgets and differs substantially from this supported-target summary.

## Metric availability

{table(availability)}

Undefined metrics stay undefined: reach requires a nonempty real DE budget, and the expression ratio requires a resolved positive corrected denominator. Pandas aggregation averages defined panels/strata/families separately for each metric. Thus a nominal equal-family proxy may include fewer than three families when an entire family's utility is undefined. The table and raw panel CSVs expose those omissions; summary means must not be interpreted as complete three-lineage validation or official anchors. Calibration and confirmation are different targets and cannot be compared as the same test.

## Neutral checks

{table(neutral)}

## Coverage

{table(pd.DataFrame(coverage))}

Public raw release selection contains sampled HCT/HEK treated cells, 25 H1 source targets and incomplete K562 gene coverage. Existing source provenance records quantify these limitations. Technical strata control batches; they do not identify hidden context lineage, timing, knockdown efficiency or platform. A confirmed public proxy improvement would still require a leaderboard evaluation. A failed confirmation is preserved openly and cannot be interpreted as demonstrated improvement.

## Reproduce

Use the existing pinned source caches, reviewed upstream clone and completed V8 frozen plan/model. `run_v8_pipeline.py --resume` reconstructs its prerequisites before V9. Root scripts retain every step. To reproduce from scratch, preserve the frozen experiment directory under a new name first; scripts refuse to overwrite frozen decisions. Then:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe run_v9_pipeline.py --workers 4
```

`--resume` preserves completed stages. The pipeline runs scientific tests, reserves fresh targets, evaluates calibration, tests disjoint challenge NTCs, freezes selection, evaluates confirmation once, freezes production dependencies, audits splits, preflights 24 groups, generates 900 groups, runs official prep and independently checks the full payload. It does not upload.

The independent audit checks all 360,000 source IDs, guide/depth records and integer count libraries; compares every packaged sparse value/index/pointer and axis against the generated matrix; and replays 12 complete groups (88,958,400 gene counts) exactly. Model and package SHA-256 plus code dependency hashes are saved in `analysis/`. `audit_v9_provenance.py` separately verifies linked calibration/confirmation/neutral record hashes, source file hashes, prior model references and all V7/V8 target exclusions without changing frozen numerical inputs. `.gitattributes` preserves the exact code and experiment-record bytes rather than rewriting line endings.

## Submit from Windows Command Prompt

```bat
cd /d L:\promec\Animesh\Download\Virtual
{command}
```

If interrupted after upload, inspect the returned entry ID with `vcc status ENTRY_ID --wait`; do not create a duplicate upload while the scoring job is active.
"""
    (ROOT/'V9_SUBMISSION.md').write_text(text,encoding='utf-8')
    readme=ROOT/'README.md';old=readme.read_text(encoding='utf-8');marker='<!-- V9_VERIFIED_SUBMISSION -->'
    if marker not in old:
        block='\n'+marker+f'\n\n## V9 verified upload candidate\n\n`experiments/v9_activation_preserving/prediction.vcc` is ready: 360,000 cells x 18,533 genes, all 900 target/context groups. Official prep and independent full-payload audit passed. Frozen decoder `{model["generator"]}` initializes from simulated perturbations; V4 biological response effects remain unchanged. Public confirmation improvement: **{result["confirmed_improvement"]}** (utility difference {result["utility_difference"]:+.6f}); no hidden score is claimed. V4 remains the best published submission. See [V9_SUBMISSION.md](V9_SUBMISSION.md) for every metric, limitations, reproducible steps and the CMD upload command. V8 is retained as a documented full-moment pilot.\n'
        old=old.replace('## Completed work',block+'\n\n## Completed work',1)
        readme.write_text(old,encoding='utf-8')
    if (A/'leaderboard_result.json').exists():
        from record_v9_result import document_publication
        document_publication()
    print('V9 documentation written after verified packaging',flush=True)

if __name__=='__main__':main()
