"""Record published V10 with read-only status and preserve all frozen generation evidence."""
import argparse,json,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import build_submission as build
from record_v5_result import RAW
ROOT=build.ROOT;A=ROOT/'experiments/v10_population_transport/analysis'
ENTRY='CYpAqJP7bUYqjA2ouD7U'
REPORTED={'overall':.11675412136452572,'pds':.5724641018827651,'mse':.049814892860030824,
 'nmae':.04988569676208535,'fid':-.022334971298256,'reach':.04654672162044448,'jac':.004148286360084445}
SNAPSHOTS=['submission_ready.json','independent_validation.json','reference_model.json','count_selection.json',
 'confirmation_result.json','validation_plan.json','decoder_neutral_manifest.json','linked_provenance_validation.json',
 'corrected_depth_validation.json','production_depth_audit.json','memory_safe_stage.json','runtime_versions.json']


def table(frame):
    columns=list(frame.columns);lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
    for _,row in frame.iterrows():
        values=[]
        for value in row:
            values.append('undefined' if pd.isna(value) else f'{value:.6f}' if isinstance(value,float) else str(value))
        lines.append('| '+' | '.join(values)+' |')
    return '\n'.join(lines)+'\n'


def overlay(text,block):
    begin,end='<!-- V10 publication begin -->','<!-- V10 publication end -->'
    if begin in text:
        start=text.index(begin);stop=text.index(end,start)+len(end)
        text=text[:start]+text[stop:].lstrip('\n')
    title,rest=text.split('\n',1)
    return title+'\n\n'+begin+'\n\n'+block+'\n\n'+end+'\n\n'+rest.lstrip('\n')


def document_publication():
    result=json.loads((A/'leaderboard_result.json').read_text(encoding='utf-8'))
    scaled=pd.read_csv(A/'leaderboard_comparison.csv');raw=pd.read_csv(A/'raw_metric_comparison.csv')
    block=(f"Published as `{ENTRY}`, overall **{result['overall']:.9f}**, rank **{result['rank_at_user_report']}** in the user report "
           f"and **{result['rank_at_api_check']}** at the recorded API check. V10 is the best published project submission, "
           f"improving on V4 by **{result['overall_change_from_v4']:+.9f}**. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md).")
    path=ROOT/'V10_SUBMISSION.md';text=overlay(path.read_text(encoding='utf-8'),block)
    text=text.replace('This is an unsubmitted experimental prediction. V4 remains the best published submission until another official result beats it.',
                      'This experimental prediction is now published and is the best recorded project submission. Public development tests below remain separate from official leaderboard evidence.')
    text=text.replace('- Uploaded: false. Hidden leaderboard performance: unknown.',
                      '- Published: true. Generation-time upload=false records are preserved; the publication overlay records the official score.')
    text=text.replace('## Submit from Windows Command Prompt','## Historical upload command (already published)')
    text=text.replace('If Ctrl+C interrupts waiting after upload, use `vcc status <entry-id> --wait` to follow the existing scoring job instead of creating a duplicate.',
                      f'This entry is already published. Inspect it with `vcc status {ENTRY}`. Repeating the historical submit command creates a duplicate entry.')
    path.write_text(text.rstrip()+'\n',encoding='utf-8')
    report=f'''# V10 published leaderboard result

{block.replace(' See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md).', '')} Ranks change over time. V4 and V10 use the same partition, panel and anchors.

## Scaled score comparison

{table(scaled[['metric','v4','v9','v10','change_from_v4','overall_contribution_change_from_v4']])}
Four scaled components improve versus V4: discrimination, expression accuracy, direction fidelity and significance overlap. Fold-change accuracy and reach decline. Discrimination contributes **{result['pds_overall_contribution']:+.9f}** and fidelity **{result['fidelity_overall_contribution']:+.9f}** to the overall gain. Together they account for **{result['pds_fidelity_fraction_of_gain']:.2%}** of the net improvement. Contributions are arithmetic changes in the six-metric mean, not causal attribution or significance estimates. Scaled MSE/NMAE are accuracy scores: higher is better; V4's zero expression score was clipping, not zero error.

## Raw metric comparison

{table(raw[['metric','v4','v10','change_from_v4','improved_from_v4']])}
Raw fidelity yield rises **{result['raw_fidelity_relative_change_from_v4']:.2%}**, expression error falls **{-result['raw_expression_relative_change_from_v4']:.2%}**, and reach falls **{-result['raw_reach_relative_change_from_v4']:.2%}**. Direction fidelity remains below the anchor level (its scaled component is still negative). Raw metrics and anchored scores are different quantities.

## What the result establishes

V10 beats V4 on the scored validation panel. It changes both the response estimator (five-source pooling, bounded knockdown-dose correction and uncertainty shrinkage) and count emission (compact synthetic populations with two-moment integer fitting). The result does not isolate either change's contribution. Corrected public validation anticipated the average fidelity gain and the reach decline, but source-family behavior was uneven; HCT116 fidelity regressed with full references. Compact emission still has substantial neutral-control DE bias. A leaderboard improvement does not establish faithful biological single-cell variance.

Keep the exact V10 artifact and frozen parameters as the current published baseline. Next compare V4 responses with the same compact emitter to isolate response pooling from emission; reserve any amplitude/dispersion or coverage choices for a new experiment. Aim to recover reach and fold-change accuracy while protecting fidelity and discrimination. Do not select another model using only expression error, or infer tissue-specific causes from aggregate scores. One submission cannot quantify random-seed uncertainty or establish final-round generalization.

## Provenance and reproduction

Entry submitted `{result['submitted_utc']}`; user-reported upload duration 17 minutes 5 seconds. Partition `{result['partition']}`, panel `{result['panel']}`, anchors `{result['anchor_set']}`. The independently verified local package SHA-256 is `{result['package_sha256']}`. The status API does not expose the uploaded package hash; local integrity and publication are recorded separately.

Generation-time upload=false snapshots, model code, choices and validation records remain unchanged. `publication_manifest.json` links those snapshots with the saved API response, result, comparisons and publication documents. Large prediction files remain excluded from Git.

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
.\\.competition_env\\Scripts\\python.exe record_v10_result.py
```

For offline reproduction from the saved status response, add `--api-file experiments/v10_population_transport/analysis/published_api_status.json`. These commands inspect/document the existing entry; they do not submit.
'''
    (ROOT/'V10_LEADERBOARD_RESULTS.md').write_text(report.rstrip()+'\n',encoding='utf-8')
    path=ROOT/'README.md';text=path.read_text(encoding='utf-8')
    text=text.replace('Seven submissions (first, V2, V4, V5, V6, V7 and V9) have been published.',
                      'Eight submissions (first, V2, V4, V5, V6, V7, V9 and V10) have been published.')
    text=text.replace('V4 is the best so far on the validation panel, with overall score -0.050674; direction fidelity and expression error remain weaknesses.',
                      'V10 is the best recorded validation submission, with overall score **0.116754121**, improving on V4 by **0.167428032**. Reach and fold-change accuracy remain weaknesses.')
    text=text.replace('V9 is now published below; V4 remains the preferred published reference.', 'V9 is published below; V10 is now the best published reference.')
    text=text.replace('V4 remains best. See [V7_LEADERBOARD_RESULTS.md]', 'V4 was best at V7 publication. See [V7_LEADERBOARD_RESULTS.md]')
    text=text.replace('Keep V4 as best. See [V9_LEADERBOARD_RESULTS.md]', 'V4 remained best at V9 publication; V10 is now best. See [V9_LEADERBOARD_RESULTS.md]')
    text=text.replace('See [V4_SUBMISSION.md](V4_SUBMISSION.md) for the published comparison and reproducibility.', 'See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md) for the current result and [V10_SUBMISSION.md](V10_SUBMISSION.md) for reproducibility.')
    text=text.replace('V4 remains best at **-0.050674**.', 'V4 was best at that publication, at **-0.050674**.')
    text=text.replace('V4 remains best.', 'V4 was best at that publication.')
    text=text.replace('V4 remains the best of six published project submissions.', 'At that stage, V4 was the best of six published project submissions.')
    text=text.replace('V4 is the best of the six published project submissions on this validation panel;', 'At that stage, V4 was the best of the six published project submissions on this validation panel;')
    text=text.replace('Overall remains below the zero-score anchor, with direction fidelity the largest negative component.', 'At V4 publication, overall was below the zero-score anchor, with direction fidelity the largest negative component.')
    text=text.replace('V4 remains the best published project candidate.', 'V4 was the best published candidate when this October 7 plan was written.')
    text=text.replace('V4 remains the starting point.', 'That earlier V6 strategy used V4 as its starting point.')
    rows=text.splitlines()
    for i,row in enumerate(rows):
        if row.startswith('The latest published submission is'):
            rows[i]=f"The latest and best published submission is V10 (`{ENTRY}`), overall **{result['overall']:.9f}**. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md)."
        if row=='## V10 upload candidate':rows[i]='## V10 published result'
        if row.startswith('V10 is ready at') or row.startswith('V10 is published as'):
            rows[i]=f"V10 is published as `{ENTRY}`, overall **{result['overall']:.9f}**, an improvement of **{result['overall_change_from_v4']:+.9f}** over V4. Official prep and independent full-payload validation passed. The pooled dose-corrected response and compact emitter improve discrimination and fidelity, while reach and fold-change accuracy decline. Compact neutral-control bias remains. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md), [V10_STRATEGY.md](V10_STRATEGY.md) and [V10_SUBMISSION.md](V10_SUBMISSION.md)."
    path.write_text('\n'.join(rows).rstrip()+'\n',encoding='utf-8')
    for name,marker in [('PERFORMANCE_IMPROVEMENT_PLAN.md','<!-- V10 published decision -->'),('V10_STRATEGY.md','<!-- V10 published result -->')]:
        path=ROOT/name;text=path.read_text(encoding='utf-8')
        if marker not in text:
            text=text.rstrip()+f"\n\n{marker}\n\n## V10 publication update\n\nV10 is published at **{result['overall']:.9f}**, beating V4 by **{result['overall_change_from_v4']:+.9f}**, and is now the best project submission. Discrimination and fidelity account for {result['pds_fidelity_fraction_of_gain']:.2%} of the net gain; expression and overlap also improve, while reach and fold-change accuracy worsen. Preserve V10 and use controlled ablations to separate response pooling from emission. Earlier V4-best statements describe historical decisions. Neutral-control bias remains unresolved. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md).\n"
        path.write_text(text.rstrip()+'\n',encoding='utf-8')
    artifacts=['leaderboard_result.json','published_api_status.json','leaderboard_comparison.csv','raw_metric_comparison.csv']
    build.write_json(A/'publication_manifest.json',{'entry_id':ENTRY,'recorder_sha256':build.sha256(ROOT/'record_v10_result.py'),
        'publication_artifact_sha256':{name:build.sha256(A/name) for name in artifacts},
        'publication_document_sha256':build.sha256(ROOT/'V10_LEADERBOARD_RESULTS.md'),
        'submission_document_sha256':build.sha256(ROOT/'V10_SUBMISSION.md'),
        'generation_snapshot_sha256':{name:build.sha256(A/name) for name in SNAPSHOTS},
        'scope':'Publication overlay; frozen generation/calibration/confirmation and corrected diagnostics unchanged.'})


def main(api_file=None):
    before={name:build.sha256(A/name) for name in SNAPSHOTS}
    if api_file:api=json.loads(Path(api_file).read_text(encoding='utf-8-sig'))
    else:
        process=subprocess.run([sys.executable,'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45)
        api=json.loads(process.stdout)
    assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
    for key,value in REPORTED.items():np.testing.assert_allclose(api['score_avg' if key=='overall' else 'score_'+key],value,rtol=0,atol=1e-12)
    directories={'first':'analysis','v2':'experiments/v2_cd4_hybrid/analysis','v4':'experiments/v4_lineage_response/analysis',
                 'v5':'experiments/v5_atlas_calibrated/analysis','v6':'experiments/v6_fidelity_calibrated/analysis',
                 'v7':'experiments/v7_multisource_adaptive/analysis','v9':'experiments/v9_activation_preserving/analysis'}
    previous={name:json.loads((ROOT/directory/'leaderboard_result.json').read_text(encoding='utf-8')) for name,directory in directories.items()}
    for old in previous.values():assert (api['partition'],api['panel_id'],api['anchor_version'])==(old['partition'],old['panel'],old['anchor_set'])
    scaled=pd.DataFrame([{'metric':key,**{name:old['overall'] if key=='overall' else old['metrics'][key] for name,old in previous.items()},
        'v10':value,'change_from_v4':value-(previous['v4']['overall'] if key=='overall' else previous['v4']['metrics'][key]),
        'overall_contribution_change_from_v4':None if key=='overall' else (value-previous['v4']['metrics'][key])/6} for key,value in REPORTED.items()])
    delta=REPORTED['overall']-previous['v4']['overall'];fid=(REPORTED['fid']-previous['v4']['metrics']['fid'])/6;pds=(REPORTED['pds']-previous['v4']['metrics']['pds'])/6
    np.testing.assert_allclose(scaled.overall_contribution_change_from_v4.sum(),delta,rtol=0,atol=1e-12)
    np.testing.assert_allclose(np.mean([v for k,v in REPORTED.items() if k!='overall']),REPORTED['overall'],rtol=0,atol=1e-12)
    raw=pd.DataFrame([{'metric':key,'v4':previous['v4']['raw_metrics'][key],'v9':previous['v9']['raw_metrics'][key],'v10':api[key],
        'change_from_v4':api[key]-previous['v4']['raw_metrics'][key],
        'improved_from_v4':api[key]>previous['v4']['raw_metrics'][key] if higher else api[key]<previous['v4']['raw_metrics'][key]} for key,higher in RAW.items()])
    ready=json.loads((A/'submission_ready.json').read_text(encoding='utf-8'));audit=json.loads((A/'independent_validation.json').read_text(encoding='utf-8'))
    assert audit['passed'] and audit['package_sha256']==ready['sha256']
    model=json.loads((A/'reference_model.json').read_text(encoding='utf-8'));assert model['candidate']=='dose_compact'
    for name,digest in model['dependencies'].items():assert build.sha256(ROOT/name)==digest,name
    relative=lambda key:api[key]/previous['v4']['raw_metrics'][key]-1
    result={'source':'User-reported CLI publication independently confirmed with read-only vcc --json status.',
        'recorded_utc':datetime.now(timezone.utc).isoformat(),'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published',
        'submitted_utc':api['submission_date'],'rank_at_user_report':585,'rank_at_api_check':api['rank'],
        'overall':REPORTED['overall'],'metrics':{k:v for k,v in REPORTED.items() if k!='overall'},'raw_metrics':{key:api[key] for key in RAW},
        'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],
        'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':1025,
        'package_hash_scope':'Independently verified local artifact; the status API does not expose the uploaded file hash.',
        'overall_change_from_v4':delta,'fidelity_overall_contribution':fid,'pds_overall_contribution':pds,'pds_fidelity_fraction_of_gain':(fid+pds)/delta,
        'raw_fidelity_relative_change_from_v4':relative('de_wilcoxon_direction_fidelity_yield_raw'),
        'raw_expression_relative_change_from_v4':relative('expr_mse_unbiased_capped_norm'),
        'raw_reach_relative_change_from_v4':relative('de_wilcoxon_direction_reach_raw'),
        'interpretation':'V10 is best published project candidate. Four scaled metrics improve; reach and fold-change accuracy worsen. Response and emission change together; causal attribution needs ablation. Neutral-control bias remains.',
        'generation_snapshot_sha256':before,'recorder_sha256':build.sha256(ROOT/'record_v10_result.py')}
    build.write_json(A/'published_api_status.json',api);build.write_json(A/'leaderboard_result.json',result)
    scaled.to_csv(A/'leaderboard_comparison.csv',index=False);raw.to_csv(A/'raw_metric_comparison.csv',index=False)
    document_publication();assert before=={name:build.sha256(A/name) for name in SNAPSHOTS}
    print(scaled[['metric','v4','v10','change_from_v4','overall_contribution_change_from_v4']].to_string(index=False))
    print('V10 publication documented; all frozen generation and validation records unchanged.',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--api-file',type=Path)
    main(parser.parse_args().api_file)
