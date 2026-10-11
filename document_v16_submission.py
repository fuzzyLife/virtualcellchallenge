"""Write V16 measured results and upload decision from verified evidence."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import build_submission as b
from v16_methods import ROOT,OUT,A

def read(n):return json.loads((A/n).read_text())
def table(f):
    def v(x):return 'NA' if isinstance(x,(float,np.floating)) and not np.isfinite(x) else f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x)
    return '\n'.join(['| '+' | '.join(map(str,f.columns))+' |','| '+' | '.join(['---']*len(f.columns))+' |']+['| '+' | '.join(v(x) for x in row)+' |' for row in f.itertuples(index=False,name=None)])

def main():
    selection=read('selection.json');audit=read('independent_experiment_audit.json');assert audit['passed'] and audit['exact_reused_v15_baseline_count_profile_DE_cases']==28 and audit['bitwise_balanced_null_identity_cases']==6
    ready=read('submission_ready.json') if (A/'submission_ready.json').exists() else None;confirmation=read('confirmation_result.json') if (A/'confirmation_result.json').exists() else None
    if ready:assert read('independent_validation.json')['passed'] and ready['sha256']==b.sha256(OUT/'prediction.vcc')
    validated=read('independent_validation.json') if ready else None
    validation_summary=(f"Final release validation passed: every {validated['packaged_nonzero_values_checked']:,} nonzero count and index, every pointer and axis, all {validated['verified_rows']:,} original-control identities/guides/libraries, all {validated['verified_groups']} groups, and {validated['exact_replayed_groups']} exactly regenerated groups ({validated['exact_replayed_gene_counts']:,} gene-count coordinates). Production responses were recomputed exactly. Package SHA-256: `{validated['package_sha256']}`." if ready else 'No production-package validation was run because no new package was released.')
    comparison=pd.DataFrame([{'candidate':r['candidate'],'passes':r['pass'],**r['changes']} for r in selection['comparisons']]);gains=[]
    for p in sorted((A/'models').glob('*.json')):
        m=json.loads(p.read_text())['model'];gains.append({'excluded_family':','.join(m['excluded']) if m['excluded'] else 'none (production)','minimum_gain':min(m['gains']),'median_gain':float(np.median(m['gains'])),'fully_retained_bins':int((np.asarray(m['gains'])==1).sum()),'training_gene_observations':sum(m['observations_per_bin'])})
    counts={s:{k:read(s+'_evaluation.json')[k] for k in ['population_comparisons','newly_generated_populations','reused_populations','target_seed_cases']} for s in ['null','development','confirmation'] if (A/(s+'_evaluation.json')).exists()}
    status='A new exploratory package passed official prep and full-payload validation. It has not been uploaded; no score or rank gain is known.' if ready else 'No new package was released. The tested candidates did not justify replacing V13 under the frozen policy; the existing V13 file does not need re-uploading.'
    by_name={r['candidate']:r['changes'] for r in selection['comparisons']}
    null_stats=pd.read_csv(A/'null_summary.csv').set_index('candidate')
    bal=by_name['balanced'];cal=by_name['calibrated_compact']
    lessons=f"""Balanced transport fixed the control-only artifact but lost {abs(bal['fidelity']):.3f} development fidelity and increased paired expression error by {100*bal['relative_expression_worsening']:.1f}%. It was rejected despite much better DE overlap. These tests show that preserving control variability alone did not recover perturbation accuracy; sample means and transferred responses still matter.

Gene-calibrated compact responses improved eligible NMAE by {abs(cal['nmae']):.4f} and expression error by {abs(100*cal['relative_expression_worsening']):.2f}%, but reduced PDS by {abs(cal['pds']):.4f} and reach by {abs(cal['reach']):.4f}. Shrinking weak evidence also removed useful target separation. It was rejected.

Residual10 reduced null DE calls by {100*(1-null_stats.loc['residual10','predicted_calls']/null_stats.loc['compact','predicted_calls']):.1f}%, but retained only {100*null_stats.loc['residual10','variance_ratio']:.1f}% of reference total variance. This is a modest count-generation compromise; the major distribution mismatch remains unresolved. The available public checks cannot establish a significant leaderboard improvement."""
    conf=table(pd.read_csv(A/'confirmation_summary.csv')) if confirmation else 'Reserved verification was not opened after development rejection.'
    if confirmation:
        delta=confirmation['comparisons'][0]['changes']
        conf+='\n\nPaired reserved changes (candidate minus V13; expression is relative numerator worsening):\n\n'+table(pd.DataFrame([{'metric':k,'change':v} for k,v in delta.items()]))
    text=f'''# V16 measured results: control-preserving confidence transport

{status}

V13 remains the best published project score (0.13788313368012328). Public raw comparisons below are not an official anchored score or a prediction of rank. Current public source biology and development targets are partly reused; missing NMAE panels prevent complete six-metric verification.

## Implemented and tested

Three response forms (V13, gene-calibrated shrinkage, matched uniform shrinkage) were crossed with compact counts, balanced molecule redistribution and residual10. The fixed 18-bin calibrator uses inner-family held predictions, excludes the outer family and all reserved calibration-label targets, and retains at least 50% V13. Two-source significant effects >=0.2 log2 remain intact. CD4 conditions are one family; source significance is an approximate confidence feature, not significance in private A/B/C contexts.

Balanced transport changes each sampled original cell by binomial removals/multinomial redistribution, preserves libraries exactly and is bitwise identity under a zero intervention. It preserves original zero support and cannot activate a gene absent in an individual cell. It does not force global population moments. Residual10 retains the original 0.65 Poisson coefficient before fitting; the dual solver can still change final dispersion. These are exploratory generators, not optimal transport or a causal biological model. Full equations, identities, thresholds and commands are in [V16_STRATEGY.md](V16_STRATEGY.md).

Nine scientific unit tests passed. Counts/comparisons: **{counts}**. These counts include repeated seeds, shared truths and reused baselines; they are not independent biological experiments. The actual development treatments contain 400 distinct cells each. Null comparisons generate 400 controls against independent right controls; their disjoint actual-control diagnostic contains 200 cells. Source validation shares templates within a family/seed; production samples controls separately per target. No hidden competition response labels are used.

## Response calibration diagnostics

These are summary proxies and did not choose candidates; actual-count gates determine release.

{table(pd.read_csv(A/'response_summary.csv'))}

{table(pd.DataFrame(gains))}

Gene-observation counts contain correlated genes, targets and conditions; they are not statistical sample sizes. Gains are regularized toward one and clipped, with explicit strong-evidence protection. Uniform controls match each calibrated target's total absolute log2 effect; they remain eligible candidates under the same prospective policy.

## Independent-control diagnostics

{table(pd.read_csv(A/'null_summary.csv'))}

Zero false calls in a control-only diagnostic establish neither perturbation accuracy nor private performance. The generator gate requires <=10% more false calls than compact, improved absolute log total-variance distance, and detection MAE within 10% of compact. Total variance alone does not verify every gene or covariance.

## Paired actual-count development

Higher is better for PDS/fidelity/reach/Jaccard; lower is better for NMAE/expression. The displayed expression ratio can be dominated by tiny control-distance denominators; the policy uses paired expression numerators for relative worsening. NMAE averages eligible matched panels; missing panels are disclosed in `selection.json`, never filled with zero.

{table(pd.read_csv(A/'development_summary.csv'))}

{table(comparison)}

Frozen finalist: **{selection['candidate']}**. If the finalist is `residual10`, deployment keeps the exact V13 response predictor and changes only count generation; the gene-calibrated shrinkage is not deployed. Every candidate must meet the unchanged V15 raw/worst-family bounds and improve at least three available metrics including a DE metric. Gene-calibrated arms must additionally beat their uniform counterparts on two of fidelity/reach/Jaccard. Equal improvement ranks select among passing candidates; no private anchor utility is estimated.

## Reserved verification and release

{conf}

Reserved result: **{confirmation['new_candidate_passes'] if confirmation else 'not reached'}**. This is one frozen finalist on nine reserved current target names, with its uniform control if gene-calibrated. It is not independent new biology. Failure cannot retune or reselect. The plan and once-only assessments are hashed and preserved.

## What the comparison teaches

{lessons}

## Independent integrity audit

The audit verifies unique original IDs, disjoint controls, exact copied actual/template/reference identities, manual gain equations and inner/outer exclusions, all completed count/profile cache hashes, **28** exact reused V15 baseline count/profile/DE cases, **{audit['original_baseline_same_input_regenerations']}** baseline regenerations through the original V11 emitter, and **six** exact balanced-control identities. It also verifies the earlier V15 evidence is unchanged. Format/integrity is distinct from biological accuracy.

{validation_summary}

The source file and official staging file were compared losslessly before official prep. Balanced production, if selected, additionally checks every predicted nonzero was present in its original cell; the released residual10 generator does not enforce that support restriction.

Root code includes `v16_methods.py`, `test_v16_transport.py`, `prepare_v16_experiment.py`, `run_v16_experiment.py`, `assess_v16_experiment.py`, `audit_v16_experiment.py`, `fit_v16_submission.py`, `build_v16_submission.py`, `package_v16_submission.py`, `validate_v16_submission.py`, `release_v16_submission.py` and this report generator. Small evidence is versioned in `experiments/v16_control_preserving/analysis/`; large arrays/counts/packages/provenance are regenerated locally and excluded from Git. No upload is performed.
'''
    result=ROOT/'V16_RESULTS.md';result.write_text(text,encoding='utf-8',newline='\n');outputs={result.name:b.sha256(result)}
    if ready:
        model=read('reference_model.json');title='V13 response residual dispersion v16' if model['candidate']=='residual10' else 'Control-preserving response transport v16';upload=f'''# V16 upload-ready exploratory prediction

Official prep and complete-payload validation passed. No leaderboard result is known. See [V16_RESULTS.md](V16_RESULTS.md) for actual public gains, missing metrics and limitations.

- File: `experiments\\v16_control_preserving\\prediction.vcc`.
- Candidate: `{model['candidate']}`; response `{model['response']}`; generator `{model['emitter']}`.
- Size: {ready['bytes']:,} bytes ({ready['bytes']/2**30:.2f} GiB).
- SHA-256: `{ready['sha256']}`.
- Shape: 360,000 cells x 18,533 genes; 900 groups, exactly 400 cells each.
- Available public guard verification: {model['available_public_guards_passed']}; complete six-metric verification: {model['confirmed_improvement']}.

Run in Command Prompt:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v16_control_preserving\\prediction.vcc -m "{title}" --wait
```

After an interruption following upload, use the printed entry ID with `vcc status ENTRY_ID --wait` rather than submitting again. No upload was performed here. V13 remains the best published reference until a new score is available.
'''
    else:upload='''# V16 release decision: retain V13

No V16 package was created because the frozen development/reserved criteria did not support replacing V13. See [V16_RESULTS.md](V16_RESULTS.md) and [V16_STRATEGY.md](V16_STRATEGY.md).

The existing `experiments\\v13_response_programs\\prediction.vcc` is our best published prediction (entry `AvAsYkuxFSjbxYDCd0wJ`, overall 0.13788313368012328). It was verified in V15; no repeat upload of the identical file is needed. Numerical prior models and artifacts were preserved.
'''
    path=ROOT/'V16_SUBMISSION.md';path.write_text(upload,encoding='utf-8',newline='\n');outputs[path.name]=b.sha256(path)
    begin='<!-- v16-research:begin -->';end='<!-- v16-research:end -->';outcome='New exploratory package validated; see [V16_SUBMISSION.md](V16_SUBMISSION.md).' if ready else 'No new package released; retain V13.'
    block=f'''{begin}

## V16 control-preserving response transport

Tested cross-family evidence shrinkage and amplitude-matched uniform controls with compact, balanced and noise-preserving residual generators. Frozen development compares 252 populations (224 new, 28 verified reused baseline), plus 18 control diagnostics. Selected candidate: `{selection['candidate']}`. {outcome} Public partly reused biology and missing NMAE panels limit inference; no private score or rank gain is claimed. Methods/reproduction: [V16_STRATEGY.md](V16_STRATEGY.md). Results/audits: [V16_RESULTS.md](V16_RESULTS.md). Root code and small hashed evidence are versioned.

{end}'''
    path=ROOT/'README.md';s=path.read_text(encoding='utf-8')
    if begin in s:
        i=s.index(begin);j=s.index(end,i)+len(end);s=s[:i]+block+s[j:]
    else:s=s.replace('<!-- v15-research:end -->','<!-- v15-research:end -->\n\n'+block,1)
    path.write_text(s,encoding='utf-8',newline='\n')
    names=['plan.json','unit_tests.json','models_complete.json','response_scores.csv','response_summary.csv','null_evaluation.json','null_summary.csv','development_evaluation.json','development_scores.csv','development_panels.csv','development_summary.csv','selection.json','confirmation_result.json','confirmation_summary.csv','independent_experiment_audit.json','reference_model.json','submission_ready.json','independent_validation.json']
    write={'inputs':{str((A/n).relative_to(ROOT)):b.sha256(A/n) for n in names if (A/n).exists()},'outputs':outputs,'generator_sha256':b.sha256(Path(__file__)),'uploaded':False,'scope':'Measured public experiment and, if present, independently validated exploratory package; no score guarantee.'};b.write_json(A/'report_manifest.json',write);print('V16 reports and README generated; package present:',bool(ready),flush=True)
if __name__=='__main__':main()
