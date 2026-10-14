"""Record published V20 metrics while preserving generation-time evidence."""
from pathlib import Path
import argparse,csv,hashlib,json,math
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis';ENTRY='B92qai3cUvEHeobzp8Fk'
SCORES={'overall':.13901156696775496,'pds':.6392880774753226,'mse':.04782300331502387,'nmae':.07957732709972225,'fid':-.016123050376923142,'reach':.06893529267065236,'jac':.014568751622731952}
RAW={'pds_cosine':True,'expr_mse_unbiased_capped_norm':False,'de_wilcoxon_lfc_nmae':False,'de_wilcoxon_direction_fidelity_yield_raw':True,'de_wilcoxon_direction_reach_raw':True,'de_wilcoxon_sig_jaccard':True}
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def close(a,b):assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12),(a,b)
def save_csv(p,rows):
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def table(rows):
 keys=list(rows[0]);lines=['| '+' | '.join(keys)+' |','| '+' | '.join(['---']*len(keys))+' |']
 for r in rows:lines.append('| '+' | '.join('' if v is None else f'{v:.6f}' if isinstance(v,float) else str(v) for v in r.values())+' |')
 return '\n'.join(lines)+'\n'
def main(api_file):
 api=read(api_file);baseline=read(ROOT/'experiments/v16_control_preserving/analysis/leaderboard_result.json');assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
 for key,value in SCORES.items():close(value,api['score_avg' if key=='overall' else 'score_'+key])
 assert (api['partition'],api['panel_id'],api['anchor_version'])==(baseline['partition'],baseline['panel'],baseline['anchor_set'])
 close(sum(SCORES[k] for k in SCORES if k!='overall')/6,SCORES['overall']);rows=[]
 for key,value in SCORES.items():
  old=baseline['overall'] if key=='overall' else baseline['metrics'][key];rows.append({'metric':key,'v16':old,'v20':value,'change_from_v16':value-old,'overall_contribution_change':None if key=='overall' else (value-old)/6})
 delta=SCORES['overall']-baseline['overall'];close(sum(r['overall_contribution_change'] or 0 for r in rows),delta);assert delta<0
 raw=[{'metric':key,'v16':baseline['raw_metrics'][key],'v20':api[key],'change_from_v16':api[key]-baseline['raw_metrics'][key],'relative_change':api[key]/baseline['raw_metrics'][key]-1,'higher_is_better':higher,'improved':api[key]>baseline['raw_metrics'][key] if higher else api[key]<baseline['raw_metrics'][key]} for key,higher in RAW.items()]
 ready=read(A/'submission_ready.json');audit=read(A/'independent_validation.json');model=read(A/'reference_model.json');assert ready['upload_ready'] and audit['passed'] and ready['sha256']==audit['package_sha256'];assert (A.parent/'prediction.vcc').stat().st_size==ready['bytes'];assert sha(A/'reference_model.json')==audit['model_sha256']
 for name,digest in model['dependencies'].items():assert sha(ROOT/name)==digest,'Changed frozen model dependency: '+name
 archives={'V20_RESULTS.md':'generation_results_document.md','V20_SUBMISSION.md':'generation_submission_document.md'}
 previous=A/'publication_manifest.json'
 if previous.exists():
  saved=read(previous);assert sha(A/'submission_ready.json')==saved['generation_receipt_sha256']
  for name,archive in archives.items():assert sha(A/archive)==saved['archived_generation_documents_sha256'][name]
 for name,archive in archives.items():
  p=A/archive
  if not p.exists():p.write_bytes((ROOT/name).read_bytes())
 result={'source':'User-reported CLI result confirmed by read-only official vcc --json status.','recorded_utc':datetime.now(timezone.utc).isoformat(),'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published','submitted_utc':api['submission_date'],'rank_at_user_report':552,'rank_at_api_check':api['rank'],'overall':SCORES['overall'],'metrics':{k:v for k,v in SCORES.items() if k!='overall'},'raw_metrics':{k:api[k] for k in RAW},'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':986,'overall_change_from_v16':delta,'relative_overall_change_from_v16':delta/baseline['overall'],'best_recorded_project_validation_score':False,'retained_best_entry_id':baseline['entry_id'],'generation_ready_sha256':sha(A/'submission_ready.json'),'recorder_sha256':sha(Path(__file__)),'interpretation':'Two scaled components improve (expression/fidelity); four decline. Reach accounts arithmetically for most of the net score decline. Public reserved improvement did not generalize to anonymous contexts; aggregate scores do not identify a per-gene causal mechanism.'}
 write(A/'leaderboard_result.json',result);save_csv(A/'leaderboard_comparison.csv',rows);save_csv(A/'raw_metric_comparison.csv',raw)
 reach=next(r for r in rows if r['metric']=='reach');rawreach=next(r for r in raw if r['metric']=='de_wilcoxon_direction_reach_raw');diagnostic=read(A/'published_effect_magnitude_diagnostic.json')
 report=f'''# V20 published leaderboard result

V20 is published as `{ENTRY}`, overall **{SCORES['overall']:.9f}**, a change of **{delta:+.9f}** ({delta/baseline['overall']:.2%}) from V16 on the same panel and anchors. **V16 remains the best recorded submission.** Rank was **552** in the user report and **{api['rank']}** at the API check; historical ranks cannot establish improvement when the leaderboard changes.

## Scaled comparison

{table(rows)}
All six displayed components are scaled scores: **higher is better**, including MSE and NMAE. Expression accuracy and direction fidelity improve; PDS, fold-change accuracy, reach and significance overlap decline. Reach contributes **{reach['overall_contribution_change']:+.9f}** to the overall change, **{reach['overall_contribution_change']/delta:.1%}** of the net decrease. This is arithmetic decomposition, not causal attribution.

## Raw comparison

{table(raw)}
Raw DE reach falls **{-rawreach['relative_change']:.2%}**. Raw expression error decreases slightly and raw fidelity improves slightly; raw NMAE worsens. The public reserve predicted reach improving by +.022051, but the hosted raw reach decreases by {rawreach['change_from_v16']:.6f}. Passing the available public guards did not establish transfer to anonymous A/B/C. Missing HCT116 NMAE and unaudited pretraining overlap remain limitations; they must not be retroactively described as complete or independent biological validation.

## Interpretation and next test

The result supports retaining V16. It does not support another blind increase in ensemble weight. A separate pre-emission audit finds total absolute eligible LFC mass changes of {diagnostic['contexts'][0]['relative_LFC_mass_change']:+.3%} in A, {diagnostic['contexts'][1]['relative_LFC_mass_change']:+.3%} in B and {diagnostic['contexts'][2]['relative_LFC_mass_change']:+.3%} in C. Thus the mean correction did not uniformly shrink response amplitude. These sums include non-DE genes and provide no hidden truth: they cannot tell which genuine responses were lost. Gene-specific calibration, source/context mismatch and nonlinear count/significance effects remain hypotheses, rather than identified causes.

Before another submission, prioritize a fresh context/assay or donor holdout, verify pretrained-study overlap where possible, and diagnose response gains/losses per context and target on accessible treatment truth. Keep V16 generation fixed while testing a narrowly defined correction against matched magnitude/shuffled controls. Use a new reserved set for any new selection; the V20 reserve has been spent. No new candidate is built by this publication recorder.

## Provenance

The official API confirms publication, score, panel and anchors. It does not provide an uploaded-file checksum; the user's CLI path ties the audited local package to the entry. Original models, predictions, validation receipts and the generation-time `submission_ready.json` remain unchanged; that sealed receipt still records uploaded=false at generation. Current publication state is in `leaderboard_result.json`. Original report documents are archived byte-for-byte and mapped by `publication_manifest.json`. After rerunning the generation report writer, rerun this recorder to restore the publication overlay.

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY} --wait
.\\.competition_env\\Scripts\\python.exe record_v20_result.py
```
'''
 (ROOT/'V20_LEADERBOARD_RESULTS.md').write_text(report,encoding='utf-8',newline='\n');base=(A/archives['V20_RESULTS.md']).read_text();(ROOT/'V20_RESULTS.md').write_text(base.rstrip()+f'\n\n## Published result\n\nPublished `{ENTRY}`, overall **{SCORES["overall"]:.9f}**, below V16 by **{-delta:.9f}** ({-delta/baseline["overall"]:.2%}). Expression and fidelity improve; reach and three other scaled components decline. V16 remains the best model. See [V20_LEADERBOARD_RESULTS.md](V20_LEADERBOARD_RESULTS.md) for exact comparisons and limitations. The preceding generation report is archived unchanged.\n',encoding='utf-8',newline='\n')
 base=(A/archives['V20_SUBMISSION.md']).read_text();base=base.replace('# V20 submission: independently validated and ready for manual upload','# V20 submission: published');start=base.index('The file has **not been uploaded**');end=base.index('V20 is the bounded',start);base=base[:start]+f'Published entry `{ENTRY}`, overall **{SCORES["overall"]:.9f}**. V16 remains the best model; see [published comparison](V20_LEADERBOARD_RESULTS.md). '+base[end:];start=base.index('## Upload from Windows CMD');end=base.index('## What passed',start);base=base[:start]+f'## Published status\n\n```bat\ncd /d L:\\promec\\Animesh\\Download\\Virtual\n.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY} --wait\n```\n\nThe audited artifact is already published. Generation-time release receipts remain immutable; current status is recorded separately in `analysis/leaderboard_result.json`.\n\n'+base[end:];(ROOT/'V20_SUBMISSION.md').write_text(base,encoding='utf-8',newline='\n')
 p=ROOT/'README.md';text=p.read_text();text=text.replace('Eleven submissions (first, V2, V4, V5, V6, V7, V9, V10, V11, V13 and V16) have been published.','Twelve submissions (first, V2, V4, V5, V6, V7, V9, V10, V11, V13, V16 and V20) have been published.');text=text.replace('Its full production package is independently validated and ready for manual upload: [V20 submission and command](V20_SUBMISSION.md). No new leaderboard score is claimed; V16 remains the best published submission.','V20 is published with overall **0.139011567**, below V16 by **0.002602290** (1.84%). Expression and fidelity improve, while DE reach falls; V16 remains the best published submission. See [V20 published results](V20_LEADERBOARD_RESULTS.md) and [artifact status](V20_SUBMISSION.md).');p.write_text(text,encoding='utf-8',newline='\n')
 publication=['published_api_status.json','leaderboard_result.json','leaderboard_comparison.csv','raw_metric_comparison.csv','published_effect_magnitude_diagnostic.json'];write(A/'publication_manifest.json',{'entry_id':ENTRY,'generation_receipt_sha256':sha(A/'submission_ready.json'),'generation_model_sha256':sha(A/'reference_model.json'),'archived_generation_documents_sha256':{name:sha(A/archive) for name,archive in archives.items()},'publication_artifact_sha256':{name:sha(A/name) for name in publication},'publication_document_sha256':{name:sha(ROOT/name) for name in ['V20_RESULTS.md','V20_SUBMISSION.md','V20_LEADERBOARD_RESULTS.md']},'recorder_sha256':sha(Path(__file__))});print('V20 publication recorded; V16 retained; overall change',delta,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--api-file',type=Path,default=A/'published_api_status.json');main(p.parse_args().api_file)
