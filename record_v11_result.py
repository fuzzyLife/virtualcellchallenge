"""Record published V11 through read-only status, preserving frozen generation evidence."""
import argparse,json,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import build_submission as build
from record_v5_result import RAW
from record_v10_result import table
ROOT=build.ROOT;A=ROOT/'experiments/v11_consensus_coverage/analysis';ENTRY='bmIAGveuEvFiYQTvz0Ma'
REPORTED={'overall':.1370039741117335,'pds':.6382673747455824,'mse':.024946411920956996,
 'nmae':.0810581335806226,'fid':-.015859493628032434,'reach':.08002167093476484,'jac':.013589747116506833}
SNAPSHOTS=['submission_ready.json','independent_validation.json','reference_model.json','count_selection.json',
 'confirmation_result.json','validation_plan.json','decoder_neutral_manifest.json','documentation_validation.json',
 'model_validation.json','memory_safe_stage.json','runtime_versions.json','run_parameters.json','prediction_validation.json',
 'development_evaluation.json','confirmation_evaluation.json','development_power_validation.json','confirmation_power_validation.json']


def document(result,scaled,raw):
    block=(f"Published as `{ENTRY}`, overall **{result['overall']:.9f}**, rank **529** in the user report "
           f"and **{result['rank_at_api_check']}** at the recorded API check. V11 is the best published project submission, "
           f"improving on V10 by **{result['overall_change_from_v10']:+.9f}**. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md).")
    original=A/'generation_submission_document.md'
    expected=json.loads((A/'documentation_validation.json').read_text())['submission_document_sha256']
    if not original.exists():
        assert build.sha256(ROOT/'V11_SUBMISSION.md')==expected
        original.write_bytes((ROOT/'V11_SUBMISSION.md').read_bytes())
    assert build.sha256(original)==expected
    text=original.read_text(encoding='utf-8');title,rest=text.split('\n',1)
    text=title+'\n\n<!-- V11 publication begin -->\n\n'+block+'\n\n<!-- V11 publication end -->\n\n'+rest.lstrip('\n')
    text=text.replace('V10 remains the best published project submission at 0.116754121 until another official result beats it. V11 has not been uploaded.',
                      'At generation time, V10 was the best published submission at 0.116754121. V11 is now published; the original generation-time document is archived in experiments/v11_consensus_coverage/analysis/generation_submission_document.md.')
    text=text.replace('- Uploaded: false. Leaderboard performance unknown.', '- Published: true. Generation-time upload=false JSON snapshots remain unchanged; the publication record supplies the official result.')
    text=text.replace('## Submit from Windows Command Prompt','## Historical upload command (already published)')
    text=text.replace("If waiting is interrupted after upload, follow the existing entry with `vcc status <entry-id> --wait` instead of creating a duplicate.",
                      f'Inspect the existing published entry with `vcc status {ENTRY}`. The historical submit command creates a duplicate entry.')
    (ROOT/'V11_SUBMISSION.md').write_text(text,encoding='utf-8')
    report=f"""# V11 published leaderboard result

{block.replace(' See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md).', '')} Ranks change over time. V4, V9, V10 and V11 use the same validation panel and anchors.

## Scaled component comparison

{table(scaled)}
Five of six components improve over V10. The expression-accuracy component declines. Scaled MSE and NMAE are accuracy scores: higher is better. Component contributions are changes divided by six and sum to the overall change; they do not isolate the effects of individual model changes.

## Raw metric comparison

{table(raw)}
Raw reach increases **{result['raw_relative_change']['de_wilcoxon_direction_reach_raw']:.2%}**, fold-change error decreases **{-result['raw_relative_change']['de_wilcoxon_lfc_nmae']:.2%}**, and expression error increases **{result['raw_relative_change']['expr_mse_unbiased_capped_norm']:.2%}**. Fidelity yield improves slightly but its scaled score remains below the anchor. Reach remains below V4's raw reach, despite the higher overall score.

## What we learned

V11 improves on V10 on the official validation panel. The two changes are consensus-supported recovery of shrunken source responses and increased synthetic-cell dispersion. Their combined result improves discrimination, fold-change accuracy, fidelity, reach and overlap, while worsening expression error. The public confirmation checks anticipated the average reach/fold-change gains and the expression-error tradeoff. A single submission cannot isolate which change caused each improvement or establish performance in unseen final-round contexts.

Preserve V11 as the new published baseline. The next controlled experiments should protect its direction and coverage gains while improving expression accuracy: separately ablate response rescue and dispersion, then test gene-specific uncertainty/dispersion or mean-versus-pooled response calibration. Choose any new parameters using held-family validation and a separate confirmation set, rather than tuning against this score alone. The substantial neutral-control false-DE bias and SLIRP power-check reach loss remain unresolved; a better leaderboard score does not establish realistic biological cell variance.

## Provenance and reproduction

Entry submitted `{result['submitted_utc']}`. User-reported upload duration: 16 minutes 36 seconds. Partition `{result['partition']}`, panel `{result['panel']}`, anchors `{result['anchor_set']}`. Local independently audited package SHA-256: `{result['package_sha256']}`. The status API does not expose the uploaded file hash; the local artifact and server publication are linked by the user's upload record, not an independently confirmed server-side hash.

All frozen numerical code, generation choices, confirmation and count-audit records remain unchanged. The original submission document is archived byte-for-byte as `experiments/v11_consensus_coverage/analysis/generation_submission_document.md`, matching its frozen documentation-validation hash. `publication_manifest.json` links the publication records and current documents separately.

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
.\\.competition_env\\Scripts\\python.exe record_v11_result.py --api-file experiments/v11_consensus_coverage/analysis/published_api_status.json
```

These commands inspect or document the existing entry. If the generation pipeline regenerates the pre-publication report, rerun the recorder to restore the publication overlay.
"""
    (ROOT/'V11_LEADERBOARD_RESULTS.md').write_text(report,encoding='utf-8')
    path=ROOT/'README.md';text=path.read_text(encoding='utf-8')
    text=text.replace('Eight submissions (first, V2, V4, V5, V6, V7, V9 and V10) have been published.',
                      'Nine submissions (first, V2, V4, V5, V6, V7, V9, V10 and V11) have been published.')
    text=text.replace('V10 is the best recorded validation submission, with overall score **0.116754121**, improving on V4 by **0.167428032**. Reach and fold-change accuracy remain weaknesses.',
                      'V11 is the best recorded validation submission, with overall score **0.137003974**, improving on V10 by **0.020249853**. Expression accuracy and biological dispersion remain weaknesses.')
    text=text.replace('See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md) for the current result and [V10_SUBMISSION.md](V10_SUBMISSION.md) for reproducibility.',
                      'See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md) for the current result and [V11_SUBMISSION.md](V11_SUBMISSION.md) for reproducibility.')
    text=text.replace('V10 is now the best published reference.', 'V11 is now the best published reference.').replace('V10 is now best.', 'V11 is now best.')
    rows=text.splitlines()
    for i,row in enumerate(rows):
        if row.startswith('The latest and best published submission is'):
            rows[i]=f'The latest and best published submission is V11 (`{ENTRY}`), overall **{result["overall"]:.9f}**. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md).'
        if row=='## V11 upload candidate':rows[i]='## V11 published result'
        if row.startswith('V11 is ready at') or row.startswith('V11 is published as'):
            rows[i]=f'V11 is published as `{ENTRY}`, overall **{result["overall"]:.9f}**, rank 529 at publication, improving on V10 by **{result["overall_change_from_v10"]:+.9f}**. Five of six scaled components improve; expression accuracy declines. Official preparation and independent full-payload/replay/provenance audits passed. Source consensus rescue and higher dispersion improve coverage together; neutral-control bias remains. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md), [V11_STRATEGY.md](V11_STRATEGY.md) and [V11_SUBMISSION.md](V11_SUBMISSION.md).'
    path.write_text('\n'.join(rows)+'\n',encoding='utf-8')


def main(api_file=None):
    before={name:build.sha256(A/name) for name in SNAPSHOTS}
    if api_file:api=json.loads(Path(api_file).read_text(encoding='utf-8-sig'))
    else:
        p=subprocess.run([sys.executable,'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45);api=json.loads(p.stdout)
    assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
    for k,v in REPORTED.items():np.testing.assert_allclose(api['score_avg' if k=='overall' else 'score_'+k],v,rtol=0,atol=1e-12)
    previous={name:json.loads((ROOT/f'experiments/{directory}/analysis/leaderboard_result.json').read_text()) for name,directory in
              [('v4','v4_lineage_response'),('v9','v9_activation_preserving'),('v10','v10_population_transport')]}
    for old in previous.values():assert (api['partition'],api['panel_id'],api['anchor_version'])==(old['partition'],old['panel'],old['anchor_set'])
    baseline=previous['v10'];change=REPORTED['overall']-baseline['overall']
    scaled=pd.DataFrame([{'metric':k,**{n:r['overall'] if k=='overall' else r['metrics'][k] for n,r in previous.items()},'v11':v,
       'change_from_v10':v-(baseline['overall'] if k=='overall' else baseline['metrics'][k]),
       'overall_contribution_change_from_v10':None if k=='overall' else (v-baseline['metrics'][k])/6} for k,v in REPORTED.items()])
    np.testing.assert_allclose(scaled.overall_contribution_change_from_v10.sum(),change,rtol=0,atol=1e-12)
    np.testing.assert_allclose(np.mean([v for k,v in REPORTED.items() if k!='overall']),REPORTED['overall'],rtol=0,atol=1e-12)
    assert sum(REPORTED[k]>baseline['metrics'][k] for k in REPORTED if k!='overall')==5
    raw=pd.DataFrame([{'metric':k,'v4':previous['v4']['raw_metrics'][k],'v10':baseline['raw_metrics'][k],'v11':api[k],
       'change_from_v10':api[k]-baseline['raw_metrics'][k],'improved_from_v10':api[k]>baseline['raw_metrics'][k] if higher else api[k]<baseline['raw_metrics'][k]} for k,higher in RAW.items()])
    ready=json.loads((A/'submission_ready.json').read_text());audit=json.loads((A/'independent_validation.json').read_text());model=json.loads((A/'reference_model.json').read_text())
    assert audit['passed'] and audit['package_sha256']==ready['sha256'] and model['candidate']=='consensus20_noise65'
    assert (A.parent/'prediction.vcc').stat().st_size==ready['bytes']
    for name,digest in model['dependencies'].items():assert build.sha256(ROOT/name)==digest,name
    result={'source':'User-reported CLI publication confirmed with read-only vcc --json status.','recorded_utc':datetime.now(timezone.utc).isoformat(),
       'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published','submitted_utc':api['submission_date'],
       'rank_at_user_report':529,'rank_at_api_check':api['rank'],'overall':REPORTED['overall'],'metrics':{k:v for k,v in REPORTED.items() if k!='overall'},
       'raw_metrics':{k:api[k] for k in RAW},'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],
       'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':996,'overall_change_from_v10':change,
       'raw_relative_change':{k:api[k]/baseline['raw_metrics'][k]-1 for k in RAW},'generation_snapshot_sha256':before,'recorder_sha256':build.sha256(ROOT/'record_v11_result.py'),
       'interpretation':'Best published project submission. Five components improve over V10; expression accuracy declines. Response and dispersion change together; neutral-control bias remains.'}
    for path in [ROOT/'analysis/leaderboard_result.json',*ROOT.glob('experiments/v*/analysis/leaderboard_result.json')]:
        if path.parent==A:continue
        older=json.loads(path.read_text());assert result['overall']>older['overall']
    build.write_json(A/'published_api_status.json',api);build.write_json(A/'leaderboard_result.json',result)
    scaled.to_csv(A/'leaderboard_comparison.csv',index=False);raw.to_csv(A/'raw_metric_comparison.csv',index=False)
    document(result,scaled,raw);assert before=={name:build.sha256(A/name) for name in SNAPSHOTS}
    files=['published_api_status.json','leaderboard_result.json','leaderboard_comparison.csv','raw_metric_comparison.csv','generation_submission_document.md']
    docs=['V11_LEADERBOARD_RESULTS.md','V11_SUBMISSION.md']
    build.write_json(A/'publication_manifest.json',{'entry_id':ENTRY,'recorder_sha256':build.sha256(ROOT/'record_v11_result.py'),
       'publication_artifact_sha256':{n:build.sha256(A/n) for n in files},'publication_document_sha256':{n:build.sha256(ROOT/n) for n in docs},
       'generation_snapshot_sha256':before,'original_documentation_record_points_to':'generation_submission_document.md',
       'comparison_baseline_sha256':{n:build.sha256(ROOT/f'experiments/{d}/analysis/leaderboard_result.json') for n,d in [('v4','v4_lineage_response'),('v9','v9_activation_preserving'),('v10','v10_population_transport')]},
       'recorder_dependency_sha256':{n:build.sha256(ROOT/n) for n in ['build_submission.py','record_v5_result.py','record_v10_result.py']},
       'scope':'Publication overlay only. Frozen numerical dependencies, generation, selection, confirmation and package unchanged.'})
    print(scaled[['metric','v10','v11','change_from_v10','overall_contribution_change_from_v10']].to_string(index=False))
    print('V11 publication recorded; frozen generation evidence unchanged.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--api-file',type=Path);main(p.parse_args().api_file)
