"""Record the existing published V9 entry; never submit or alter frozen generation records."""
import argparse,json,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import build_submission as build
from record_v5_result import RAW
from document_v8_submission import table
ROOT=build.ROOT;A=ROOT/'experiments/v9_activation_preserving/analysis';ENTRY='QBpef7kQNiFa4IFSxNi5'
REPORTED={'overall':-.07044728067144884,'pds':.3133011753390615,'mse':.03371872311541998,
 'nmae':.07492664275951895,'fid':-.9211825122183547,'reach':.10203419909697964,'jac':-.02548191212131842}
SNAPSHOTS=['submission_ready.json','independent_validation.json','reference_model.json','count_selection.json',
 'confirmation_result.json','validation_plan.json','decoder_neutral_manifest.json','linked_provenance_validation.json']


def overlay(text,block):
    begin,end='<!-- V9 publication begin -->','<!-- V9 publication end -->'
    if begin in text:
        start=text.index(begin);stop=text.index(end,start)+len(end);text=text[:start]+text[stop:].lstrip('\n')
    title,rest=text.split('\n',1)
    return title+'\n\n'+begin+'\n\n'+block+'\n\n'+end+'\n\n'+rest.lstrip('\n')


def document_publication():
    result=json.loads((A/'leaderboard_result.json').read_text());scaled=pd.read_csv(A/'leaderboard_comparison.csv');raw=pd.read_csv(A/'raw_metric_comparison.csv')
    block=f'Published as `{ENTRY}`, overall **{result["overall"]:.9f}**, rank {result["rank_at_user_report"]} in the user report. V4 remains best at **-0.050673910**. Five scaled metrics improved versus V4, but the fidelity loss outweighed their gains. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md).'
    path=ROOT/'V9_SUBMISSION.md';text=overlay(path.read_text(),block)
    text=text.replace('Upload-ready:', 'Verified published package:').replace('No upload or leaderboard score yet.','Published leaderboard scores are recorded above; the public-development tables below retain their original meaning.')
    text=text.replace('## Submit from Windows Command Prompt','## Historical upload command (already published)')
    text=text.replace('If interrupted after upload, inspect the returned entry ID with `vcc status ENTRY_ID --wait`; do not create a duplicate upload while the scoring job is active.',f'This upload is complete. Inspect the existing result with `vcc status {ENTRY}`. Repeating the submit command creates a duplicate submission.')
    path.write_text(text,encoding='utf-8')
    comparison=scaled[['metric','v4','v7','v9','change_from_v4','overall_contribution_change_from_v4']]
    raw_table=raw[['metric','v4','v9','change_from_v4','improved_from_v4']]
    summary=pd.read_csv(A/'count_confirmation_summary.csv').set_index('candidate')
    f0=summary.loc['legacy','fidelity'];f1=summary.loc['legacy_fit50','fidelity']
    text=fr"""# V9 published leaderboard result

Entry `{ENTRY}` is **published**, overall **{result['overall']:.9f}**, rank **{result['rank_at_user_report']}** in the user report and **{result['rank_at_api_check']}** at the read-only API check. Ranks change over time. V4 remains the best published project submission at **-0.050673910**. Both entries use the same partition, panel and anchors.

## Scaled score comparison

{table(comparison)}

Five of six scaled metrics improve versus V4. Fidelity contributes **{result['fidelity_overall_contribution']:+.9f}** to the overall difference, while the other five combined contribute **{result['other_metrics_overall_contribution']:+.9f}**. Their sum is **{result['overall_change_from_v4']:+.9f}**. These are arithmetic contributions to the six-metric mean, not causal attribution or significance estimates.

## Raw metrics

{table(raw_table)}

Raw expression error falls **{abs(result['raw_expression_relative_change_from_v4']):.2%}**; raw fidelity yield falls **{abs(result['raw_fidelity_relative_change_from_v4']):.2%}**. The positive scaled expression score indicates improved accuracy relative to the anchors; V4's zero scaled expression score was clipping, not zero error. Raw discrimination, LFC error, reach and overlap also improve. A larger negative scaled fidelity component defeats those gains.

## What we learned

V9 preserved the exact V4/V7 biological effect arrays and added 50% mean/pooled-profile correction on top of V4 zero-aware simulation, using new random draws. This is a decoder experiment, not a newly learned biological response model. The result demonstrates that improved mean expression, discrimination and neutral-control preservation do not establish perturbation direction fidelity. Aggregate scores cannot show whether the loss comes from fewer correct calls, extra incorrect calls, changed ranks/variance, or specific contexts/targets. A single upload cannot estimate seed uncertainty or isolate each emission change causally.

The public confirmation proxy improved by +0.019533 and passed its declared guard. However, public mean raw fidelity already fell from **{f0:.6f}** to **{f1:.6f}**. The guard used a reference-call>=10 budget-weighted subset, whose fidelity changed from 0.002476 to 0.002724. That very low, narrowly supported improvement did not protect the overall yield. Missing metric panels also reduced the effective support of public utility averages. The recorded confirmation acceptance remains historically true under that declared rule, but it did not establish hidden improvement.

A retrospective full mean-fidelity guard allowing at most 0.01 loss would have **rejected** V9 on public confirmation. This is an audit finding, not a preregistered rule or a retuned confirmation success. The original selection and confirmation records remain unchanged.

## Next decision

Keep V4 as the published reference. Retain corrected common-control scoring and expression jackknifes, but do not promote V9's count decoder. Before another upload, require fidelity protection both on a common supported-target set and on overall yield with explicit missingness; report real/predicted/correct DE budgets per target and equal-lineage/context results. A gain confined to the weighted subset should not override a material overall fidelity loss. Add neutral and perturbation distribution/detection checks together, and use paired seeds with reference DE budgets closer to the challenge. New thresholds/folds belong to a new experiment; V9 confirmation is now used development data.

Prioritize validated response direction and target-specific uncertainty over another unchanged-model resampling upload. Do not claim a particular context, dose or gene is responsible from these aggregate numbers. This result alone does not justify uploading a 25% decoder variant or predict that it will score better.

## Provenance and reproduce

Submitted `{result['submitted_utc']}`; upload duration 13 minutes 20 seconds from the user report. Partition `{result['partition']}`, panel `{result['panel']}`, anchors `{result['anchor_set']}`.

The independently validated local package SHA-256 is `{result['package_sha256']}`. The API does not expose the uploaded package hash; local identity and publication are recorded separately. Generation-time upload=false snapshots, code/data hashes and calibration decisions are preserved. `publication_manifest.json` records their hashes and links the publication overlay, status JSON and scaled/raw comparison CSVs.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status {ENTRY}
.\.competition_env\Scripts\python.exe record_v9_result.py
```

For an offline reproduction from the saved status JSON, add `--api-file experiments/v9_activation_preserving/analysis/published_api_status.json`. These commands inspect/document the existing entry; they do not upload. [V9_SUBMISSION.md](V9_SUBMISSION.md) retains the original numerical strategy and public evaluation results.
"""
    (ROOT/'V9_LEADERBOARD_RESULTS.md').write_text(text,encoding='utf-8')
    path=ROOT/'README.md';text=path.read_text().replace('Six submissions (first, V2, V4, V5, V6 and V7) have been published.','Seven submissions (first, V2, V4, V5, V6, V7 and V9) have been published.')
    rows=text.splitlines()
    for i,row in enumerate(rows):
        if row.startswith('The latest published submission is'):
            rows[i]=f'The latest published submission is V9 (`{ENTRY}`), overall **{result["overall"]:.9f}**. V4 remains best at **-0.050673910**. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md).'
        if row=='## V9 verified upload candidate':rows[i]='## V9 published result'
        if row.startswith('`experiments/v9_activation_preserving/prediction.vcc` is ready:') or row.startswith('V9 is published as'):
            rows[i]=f'V9 is published as `{ENTRY}`, overall **{result["overall"]:.9f}**. Official prep and independent full-payload audit passed. Its 50% count calibration improved five scaled metrics versus V4 but worsened fidelity enough to lower overall performance. V4 biological effects are unchanged. Public confirmation acceptance did not predict leaderboard improvement; its subset fidelity guard was too permissive. Keep V4 as best. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md) and [V9_SUBMISSION.md](V9_SUBMISSION.md). V8 remains an unsubmitted diagnostic pilot.'
    text='\n'.join(rows)+'\n';text=text.replace('Use the V9 candidate below for the new upload.','V9 is now published below; V4 remains the preferred published reference.')
    path.write_text(text,encoding='utf-8')
    path=ROOT/'PERFORMANCE_IMPROVEMENT_PLAN.md';text=path.read_text();marker='<!-- V9 published decision -->'
    if marker not in text:text+='\n'+marker+'\n\n## V9 publication update\n\nV9 published at -0.070447281, below V4 (-0.050673910). Five scaled metrics improved, but fidelity lost -0.237327864 and outweighed their gains. The public mean-fidelity drop was already visible; the reference-supported weighted guard was insufficient. Retain the corrected scorer, preserve V4 as best, and require both overall and supported-target fidelity protection plus explicit metric availability before promoting another decoder. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md).\n'
    path.write_text(text,encoding='utf-8')
    artifacts=['leaderboard_result.json','published_api_status.json','leaderboard_comparison.csv','raw_metric_comparison.csv']
    build.write_json(A/'publication_manifest.json',{'entry_id':ENTRY,'recorder_sha256':build.sha256(ROOT/'record_v9_result.py'),
        'publication_artifact_sha256':{n:build.sha256(A/n) for n in artifacts},
        'publication_document_sha256':build.sha256(ROOT/'V9_LEADERBOARD_RESULTS.md'),
        'submission_document_sha256':build.sha256(ROOT/'V9_SUBMISSION.md'),
        'generation_snapshot_sha256':{n:build.sha256(A/n) for n in SNAPSHOTS},
        'scope':'Publication overlay. Frozen generation and calibration/confirmation records are unchanged.'})


def main(api_file=None):
    before={n:build.sha256(A/n) for n in SNAPSHOTS}
    if api_file:api=json.loads(Path(api_file).read_text(encoding='utf-8-sig'))
    else:
        process=subprocess.run([sys.executable,'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45)
        api=json.loads(process.stdout)
    assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
    for key,value in REPORTED.items():np.testing.assert_allclose(api['score_avg' if key=='overall' else 'score_'+key],value,rtol=0,atol=1e-12)
    previous={name:json.loads((ROOT/directory/'leaderboard_result.json').read_text()) for name,directory in {
        'first':Path('analysis'),'v2':Path('experiments/v2_cd4_hybrid/analysis'),'v4':Path('experiments/v4_lineage_response/analysis'),
        'v5':Path('experiments/v5_atlas_calibrated/analysis'),'v6':Path('experiments/v6_fidelity_calibrated/analysis'),
        'v7':Path('experiments/v7_multisource_adaptive/analysis')}.items()}
    for old in previous.values():assert (api['partition'],api['panel_id'],api['anchor_version'])==(old['partition'],old['panel'],old['anchor_set'])
    scaled=pd.DataFrame([{'metric':key,**{name:old['overall'] if key=='overall' else old['metrics'][key] for name,old in previous.items()},
        'v9':value,'change_from_v4':value-(previous['v4']['overall'] if key=='overall' else previous['v4']['metrics'][key]),
        'overall_contribution_change_from_v4':None if key=='overall' else (value-previous['v4']['metrics'][key])/6} for key,value in REPORTED.items()])
    delta=REPORTED['overall']-previous['v4']['overall']
    np.testing.assert_allclose(scaled.overall_contribution_change_from_v4.sum(),delta,rtol=0,atol=1e-12)
    np.testing.assert_allclose(np.mean([v for k,v in REPORTED.items() if k!='overall']),REPORTED['overall'],rtol=0,atol=1e-12)
    raw=pd.DataFrame([{'metric':key,'v4':previous['v4']['raw_metrics'][key],'v7':previous['v7']['raw_metrics'][key],'v9':api[key],
        'change_from_v4':api[key]-previous['v4']['raw_metrics'][key],
        'improved_from_v4':api[key]>previous['v4']['raw_metrics'][key] if higher else api[key]<previous['v4']['raw_metrics'][key]} for key,higher in RAW.items()])
    ready=json.loads((A/'submission_ready.json').read_text());validation=json.loads((A/'independent_validation.json').read_text())
    assert validation['passed'] and validation['package_sha256']==ready['sha256']
    model=json.loads((A/'reference_model.json').read_text())
    assert build.sha256(A/'reference_model_arrays.npz')==model['model_arrays_sha256']
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as current, np.load(ROOT/'experiments/v7_multisource_adaptive/analysis/reference_model_arrays.npz',allow_pickle=False) as prior:
        for key in ['effects','genes','targets']:np.testing.assert_array_equal(current[key],prior[key])
    raw_expression=api['expr_mse_unbiased_capped_norm']/previous['v4']['raw_metrics']['expr_mse_unbiased_capped_norm']-1
    raw_fidelity=api['de_wilcoxon_direction_fidelity_yield_raw']/previous['v4']['raw_metrics']['de_wilcoxon_direction_fidelity_yield_raw']-1
    fid_contribution=(REPORTED['fid']-previous['v4']['metrics']['fid'])/6
    conf=pd.read_csv(A/'count_confirmation_summary.csv').set_index('candidate')
    result={'source':'User-reported CLI publication independently confirmed by read-only vcc --json status.',
        'recorded_utc':datetime.now(timezone.utc).isoformat(),'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published',
        'submitted_utc':api['submission_date'],'rank_at_user_report':1084,'rank_at_api_check':api['rank'],
        'overall':REPORTED['overall'],'metrics':{k:v for k,v in REPORTED.items() if k!='overall'},'raw_metrics':{k:api[k] for k in RAW},
        'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],
        'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':800,
        'package_hash_scope':'Independently validated local artifact; status API does not expose uploaded package hash.',
        'overall_change_from_v4':delta,'fidelity_overall_contribution':fid_contribution,'other_metrics_overall_contribution':delta-fid_contribution,
        'raw_expression_relative_change_from_v4':raw_expression,'raw_fidelity_relative_change_from_v4':raw_fidelity,
        'response_effect_arrays_equal_v4':True,'public_confirmation_original_acceptance':True,
        'retrospective_mean_fidelity_guard_passed':bool(conf.loc['legacy_fit50','fidelity']>=conf.loc['legacy','fidelity']-.01),
        'interpretation':'V4 remains best. Five scaled metric gains are outweighed by fidelity loss. The subset-weighted public guard did not adequately protect overall yield; this is a decoder experiment with unchanged biological effects.',
        'generation_snapshot_sha256':before,'recorder_sha256':build.sha256(ROOT/'record_v9_result.py')}
    build.write_json(A/'published_api_status.json',api);build.write_json(A/'leaderboard_result.json',result)
    scaled.to_csv(A/'leaderboard_comparison.csv',index=False);raw.to_csv(A/'raw_metric_comparison.csv',index=False)
    document_publication()
    assert before=={n:build.sha256(A/n) for n in SNAPSHOTS}
    print(scaled[['metric','v4','v9','change_from_v4','overall_contribution_change_from_v4']].to_string(index=False))
    print('Raw expression change:',raw_expression,'raw fidelity change:',raw_fidelity)
    print('V9 publication documented; all frozen generation/decision snapshots unchanged.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--api-file',type=Path)
    main(parser.parse_args().api_file)
