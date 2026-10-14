"""Record V24 publication without altering sealed model, package or generation receipts."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,json,math,subprocess,sys
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v24_cd4_response/analysis';ENTRY='81Y2yX7utdAau9rr3jCj'
SCORES={'overall':.14299785521745675,'pds':.644246366710111,'mse':.04756840736075477,'nmae':.08227057284911048,'fid':-.01579183745269364,'reach':.08483984010194835,'jac':.01485378173550954}
RAW={'pds_cosine':True,'expr_mse_unbiased_capped_norm':False,'de_wilcoxon_lfc_nmae':False,'de_wilcoxon_direction_fidelity_yield_raw':True,'de_wilcoxon_direction_reach_raw':True,'de_wilcoxon_sig_jaccard':True}
ARCHIVES={'V24_RESULTS.md':'generation_results_document.md','V24_SUBMISSION.md':'generation_submission_document.md'}
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def close(a,b):assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12),(a,b)
def csvfile(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def table(rows):
 keys=list(rows[0]);lines=['| '+' | '.join(keys)+' |','| '+' | '.join(['---']*len(keys))+' |']
 for r in rows:lines.append('| '+' | '.join('' if v is None else f'{v:.9f}' if isinstance(v,float) else str(v) for v in r.values())+' |')
 return '\n'.join(lines)
def main(api_file=None):
 if api_file:api=read(api_file)
 else:
  query=subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe'),'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45);api=json.loads(query.stdout);write(A/'published_api_status.json',api)
 assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
 for k,v in SCORES.items():close(v,api['score_avg' if k=='overall' else 'score_'+k])
 close(SCORES['overall'],sum(v for k,v in SCORES.items() if k!='overall')/6)
 baseline=read(ROOT/'experiments/v16_control_preserving/analysis/leaderboard_result.json');v20=read(ROOT/'experiments/v20_stack_transfer/analysis/leaderboard_result.json');assert (api['partition'],api['panel_id'],api['anchor_version'])==(baseline['partition'],baseline['panel'],baseline['anchor_set'])==(v20['partition'],v20['panel'],v20['anchor_set'])
 ready=read(A/'submission_ready.json');audit=read(A/'independent_validation.json');model=read(A/'reference_model.json');assert ready['upload_ready'] and audit['passed'] and ready['sha256']==audit['package_sha256'];assert sha(A/'independent_validation.json')==ready['independent_validation_sha256'];assert sha(A/'reference_model.json')==audit['model_sha256'] and (A.parent/'prediction.vcc').stat().st_size==ready['bytes']
 for name,h in model['dependencies'].items():assert sha(ROOT/name)==h,'Frozen model dependency changed: '+name
 preserved=['submission_ready.json','reference_model.json','independent_validation.json','independent_experiment_audit.json','plan.json','selection.json','confirmation_result.json','runtime.json'];sealed={n:sha(A/n) for n in preserved};old=A/'publication_manifest.json'
 if old.exists():
  previous=read(old);assert sealed==previous['sealed_generation_files_sha256']
  for n,p in ARCHIVES.items():assert sha(A/p)==previous['archived_generation_documents_sha256'][n]
 for n,p in ARCHIVES.items():
  if not (A/p).exists():(A/p).write_bytes((ROOT/n).read_bytes())
 rows=[]
 for k,v in SCORES.items():
  before=baseline['overall'] if k=='overall' else baseline['metrics'][k];prior=v20['overall'] if k=='overall' else v20['metrics'][k];rows.append({'metric':k,'v16':before,'v20':prior,'v24':v,'change_from_v16':v-before,'change_from_v20':v-prior,'overall_contribution_vs_v16':None if k=='overall' else (v-before)/6})
 delta=SCORES['overall']-baseline['overall'];assert delta>0 and all(r['change_from_v16']>0 for r in rows);close(sum(r['overall_contribution_vs_v16'] or 0 for r in rows),delta)
 raw=[{'metric':k,'v16':baseline['raw_metrics'][k],'v24':api[k],'change_from_v16':api[k]-baseline['raw_metrics'][k],'relative_change':api[k]/baseline['raw_metrics'][k]-1 if baseline['raw_metrics'][k] else None,'higher_is_better':higher,'improved':api[k]>baseline['raw_metrics'][k] if higher else api[k]<baseline['raw_metrics'][k]} for k,higher in RAW.items()]
 scores=[read(p)['overall'] for p in ROOT.glob('experiments/v*/analysis/leaderboard_result.json') if p.parent!=A];scores.append(read(ROOT/'analysis/leaderboard_result.json')['overall']);assert SCORES['overall']>max(scores)
 result={'source':'User-reported CLI receipt confirmed by read-only official vcc --json status. API has no uploaded-file checksum; user command links audited local artifact to entry.','recorded_utc':datetime.now(timezone.utc).isoformat(),'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published','submitted_utc':api['submission_date'],'rank_at_user_report':530,'rank_at_api_check':api['rank'],'overall':SCORES['overall'],'metrics':{k:v for k,v in SCORES.items() if k!='overall'},'raw_metrics':{k:api[k] for k in RAW},'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':978,'overall_change_from_v16':delta,'relative_overall_change_from_v16':delta/baseline['overall'],'overall_change_from_v20':SCORES['overall']-v20['overall'],'all_six_scaled_components_improve_vs_v16':True,'best_recorded_project_validation_score':True,'previous_best_entry_id':baseline['entry_id'],'generation_ready_sha256':sha(A/'submission_ready.json'),'recorder_sha256':sha(Path(__file__)),'interpretation':'Modest all-component gain from a bounded context-matched CD4 mean expert; count generator fixed and 693 unchanged groups reused. Aggregate host metrics do not identify which genes or causal mechanisms improved. Single validation-panel result does not establish statistical significance or final-test generalization.'}
 write(A/'leaderboard_result.json',result);csvfile(A/'leaderboard_comparison.csv',rows);csvfile(A/'raw_metric_comparison.csv',raw)
 summary=f'V24 is published as `{ENTRY}`, overall **{SCORES["overall"]:.9f}**, improving on V16 by **{delta:+.9f}** ({delta/baseline["overall"]:.2%}) and V20 by **{SCORES["overall"]-v20["overall"]:+.9f}**. It is the new best recorded project validation submission. All six scaled components improve versus V16.'
 report=f'''# V24 published leaderboard result

{summary}

Rank was **530** at the user report and **{api['rank']}** at the read-only API check. Historical ranks change as other teams submit; the same-panel score comparison establishes the gain.

## Scaled scores

{table(rows)}

All six displayed component scores are higher-is-better, including expression MSE and LFC NMAE. Overall is their unweighted average; the component contributions above sum exactly to the overall gain. This is arithmetic decomposition, not causal attribution. The improvement is consistent across components but modest, not a ranking breakthrough.

## Raw metrics

{table(raw)}

Raw expression MSE and raw LFC NMAE are lower-is-better; other raw metrics are higher-is-better. These values are not the CLI's scaled scores. Direction fidelity remains below its zero-score anchor; low reach and significance overlap remain substantial gaps.

## What this supports

The final expert changes 207 context-A target groups, with maximum coordinate correction0.25 log2FC, compatibility and donor-direction confidence gates, and condition weights from measured controls. The remaining693 groups reuse exactly verified V16 count blocks; original control sampling, libraries and the residual10 emitter are retained. This supports a narrow source/context calibration over another broad generic ensemble. It does not isolate the effect of individual genes, condition weights or donor gates.

Local reserve fidelity, discrimination, expression and LFC accuracy improved, while reach and overlap dipped slightly. All six hosted components now improve versus V16. The host outcome does not retroactively remove public-study overlap, missing development NMAE, shared upstream nuisance preprocessing or limited reserve-target coverage. V22 and V23 remain rejected experiments. Larger causal population-variation improvements remain unresolved.

## Provenance

Partition `{api['partition']}`, panel `{api['panel_id']}`, anchors `{api['anchor_version']}`; same as V16/V20. Submitted `{api['submission_date']}`. Upload took16 minutes18 seconds. Package is{ready['bytes']:,} bytes; SHA256 `{ready['sha256']}`. Official prep, all2,368,814,375 packaged count/index comparisons, every control identity/library check and12 exact group replays passed before submission.

The API confirms the entry and scores but exposes no uploaded-file hash. The user's upload command connects the audited local package to this entry. Generation-time uploaded=false snapshots, numerical code, model arrays and packages remain unchanged. Original report documents are archived byte-for-byte; publication_manifest.json maps generation seals and current publication evidence. Re-running the report generator restores this publication overlay from saved API evidence.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
.\\.competition_env\\Scripts\\python.exe record_v24_result.py --api-file experiments/v24_cd4_response/analysis/published_api_status.json
```
'''
 (ROOT/'V24_LEADERBOARD_RESULTS.md').write_text(report,encoding='utf-8',newline='\n');original=(A/ARCHIVES['V24_RESULTS.md']).read_text(encoding='utf-8');original=original.replace('V16 remains the best published project result, overall0.1416138568460195. V24 has not been uploaded and has no leaderboard score.','At generation time V16 was the best published score (0.1416138568460195), and V24 was not yet uploaded. The subsequent hosted result is recorded above.').replace('No automatic upload or Git push occurred. See V24_SUBMISSION.md for the exact CMD upload command.','The generator performed no automatic upload; the user subsequently submitted this exact path. Current status is in V24_SUBMISSION.md.');title,body=original.split('\n',1);(ROOT/'V24_RESULTS.md').write_text(title+'\n\n'+summary+' See [V24_LEADERBOARD_RESULTS.md](V24_LEADERBOARD_RESULTS.md).\n'+body,encoding='utf-8',newline='\n')
 submission=f'''# V24 published submission

{summary}

File: `experiments/v24_cd4_response/prediction.vcc` ({ready['bytes']:,} bytes;3.79 GiB). SHA256 `{ready['sha256']}`. Model `cdf50_de`; official packaging and independent complete-payload/identity/library/replay checks passed.

Inspect the published entry in Command Prompt:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
```

The package is already published; repeating submit would create a duplicate experiment. Original upload instructions remain in analysis/generation_submission_document.md. Methods/reproduction: V24_STRATEGY.md. Local tests: V24_RESULTS.md. Hosted comparison: V24_LEADERBOARD_RESULTS.md. Generation receipts remain immutable; current publication state is recorded separately in analysis/leaderboard_result.json.
''';(ROOT/'V24_SUBMISSION.md').write_text(submission,encoding='utf-8',newline='\n')
 f=ROOT/'README.md';text=f.read_text(encoding='utf-8-sig');start='<!-- v24-publication:begin -->';end='<!-- v24-publication:end -->';block=start+'\n\n## V24 published result\n\n'+summary+' Rank was530 at publication; the API check records'+str(api['rank'])+'. Same validation panel and anchors as V16/V20.\n\nThe bounded donor-supported CD4 expert updates207 context-A groups while693 groups retain verified V16 counts. Full package and provenance checks passed. [Published comparison](V24_LEADERBOARD_RESULTS.md), [local tests](V24_RESULTS.md), [methods](V24_STRATEGY.md), [artifact/status](V24_SUBMISSION.md). Generation-time receipts are preserved; publication is recorded separately.\n\n'+end
 if start in text:a=text.index(start);z=text.index(end,a)+len(end);text=text[:a]+block+text[z:]
 else:text=text.replace('# Virtual Cell Challenge 2026','# Virtual Cell Challenge 2026\n\n'+block,1)
 text=text.replace('Twelve submissions (first, V2, V4, V5, V6, V7, V9, V10, V11, V13, V16 and V20) have been published.','Thirteen submissions (first, V2, V4, V5, V6, V7, V9, V10, V11, V13, V16, V20 and V24) have been published.').replace('V16 is the best recorded validation submission, with overall score **0.141613857**, improving on V13 by **0.003730723** (2.71%).','V24 is the best recorded validation submission, with overall score **0.142997855**, improving on V16 by **0.001383998** (0.98%).').replace('See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md) for the current result and [V16_SUBMISSION.md](V16_SUBMISSION.md) for reproducibility.','See [V24_LEADERBOARD_RESULTS.md](V24_LEADERBOARD_RESULTS.md) for the current result and [V24_SUBMISSION.md](V24_SUBMISSION.md) for reproducibility.').replace('The new V24 file is upload-ready after official preparation and independent full-count/axis/provenance validation; it is not uploaded. V16 remains the best published submission at0.1416138568460195.','The V24 file passed official preparation and independent full-count/axis/provenance validation. It is now published at0.14299785521745675, exceeding V16 in all six scaled components.');text=text.replace('V16 remains the champion.','V16 was the champion at the time of this review.').replace('V16 remains the best published submission.','V16 remained the best published submission at that stage.').replace('V16 remains the best published reference.','V16 was the best published reference at that stage.').replace('V16 remains the best recorded submission.','V16 was the best recorded submission at that stage.').replace('It is the best recorded project validation score.','It was the best recorded project validation score at publication; V24 now exceeds it.');text=text.replace('The latest and best published submission is V16 (`wG3GGbaTsl2dXFn7cQy1`), overall **0.141613857**. See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md).','The latest and best published submission is V24 (`81Y2yX7utdAau9rr3jCj`), overall **0.142997855**. See [V24_LEADERBOARD_RESULTS.md](V24_LEADERBOARD_RESULTS.md).');f.write_text(text,encoding='utf-8',newline='\n')
 for name in ['V24_RESULTS.md','V24_SUBMISSION.md']:
  doc=ROOT/name;text_doc=doc.read_text(encoding='utf-8').replace(chr(13)+chr(10),chr(10));doc.write_bytes(text_doc.replace(chr(10),chr(13)+chr(10)).encode('utf-8'))
 assert sealed=={n:sha(A/n) for n in preserved};artifacts=['published_api_status.json','leaderboard_result.json','leaderboard_comparison.csv','raw_metric_comparison.csv'];write(A/'publication_manifest.json',{'entry_id':ENTRY,'sealed_generation_files_sha256':sealed,'archived_generation_documents_sha256':{n:sha(A/p) for n,p in ARCHIVES.items()},'publication_artifacts_sha256':{n:sha(A/n) for n in artifacts},'publication_documents_sha256':{n:sha(ROOT/n) for n in ['V24_RESULTS.md','V24_SUBMISSION.md','V24_LEADERBOARD_RESULTS.md']},'recorder_sha256':sha(Path(__file__)),'frozen_model_dependencies_verified':len(model['dependencies']),'generation_files_unchanged':True});print(table(rows));print('V24 publication recorded; all six components improve; generation seals unchanged.',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--api-file',type=Path);main(p.parse_args().api_file)
