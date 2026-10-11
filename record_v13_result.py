"""Record V13 publication without changing frozen selection or generation evidence."""
import argparse,json,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import build_submission as build
RAW={'pds_cosine':True,'expr_mse_unbiased_capped_norm':False,'de_wilcoxon_lfc_nmae':False,
 'de_wilcoxon_direction_fidelity_yield_raw':True,'de_wilcoxon_direction_reach_raw':True,'de_wilcoxon_sig_jaccard':True}

def table(frame):
 columns=list(frame.columns);lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
 for _,row in frame.iterrows():
  values=['undefined' if pd.isna(v) else f'{v:.6f}' if isinstance(v,float) else str(v) for v in row]
  lines.append('| '+' | '.join(values)+' |')
 return '\n'.join(lines)+'\n'
ROOT=build.ROOT;A=ROOT/'experiments/v13_response_programs/analysis';ENTRY='AvAsYkuxFSjbxYDCd0wJ'
REPORTED={'overall':.13788313368012328,'pds':.63913590675991,'mse':.03287835155751243,
 'nmae':.0818216118063317,'fid':-.015618423313101306,'reach':.07557363263097316,'jac':.013507722639113704}

def verify_generation(seal):
 for name,digest in seal['files'].items():
  archive=A/'generation_submission_document.md'
  path=archive if name=='V13_SUBMISSION.md' and archive.exists() else ROOT/name
  assert build.sha256(path)==digest,'Generation evidence changed: '+name
 model=json.loads((A/'reference_model.json').read_text(encoding='utf-8'))
 for name,digest in model['dependencies'].items():assert build.sha256(ROOT/name)==digest,name
 return model

def main(api_file=None):
 seal=json.loads((A/'release_manifest.json').read_text(encoding='utf-8'));seal_hash=build.sha256(A/'release_manifest.json')
 model=verify_generation(seal)
 if api_file:api=json.loads(Path(api_file).read_text(encoding='utf-8-sig'))
 else:
  query=subprocess.run([sys.executable,'-m','vcc','--json','status',ENTRY],capture_output=True,text=True,encoding='utf-8',check=True,timeout=45)
  api=json.loads(query.stdout)
 assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
 for k,v in REPORTED.items():np.testing.assert_allclose(api['score_avg' if k=='overall' else 'score_'+k],v,rtol=0,atol=1e-12)
 baseline_path=ROOT/'experiments/v11_consensus_coverage/analysis/leaderboard_result.json'
 baseline=json.loads(baseline_path.read_text(encoding='utf-8'))
 assert (api['partition'],api['panel_id'],api['anchor_version'])==(baseline['partition'],baseline['panel'],baseline['anchor_set'])
 delta=REPORTED['overall']-baseline['overall']
 scaled=pd.DataFrame([{'metric':k,'v11':baseline['overall'] if k=='overall' else baseline['metrics'][k],'v13':v,
  'change_from_v11':v-(baseline['overall'] if k=='overall' else baseline['metrics'][k]),
  'overall_contribution_change_from_v11':None if k=='overall' else (v-baseline['metrics'][k])/6} for k,v in REPORTED.items()])
 np.testing.assert_allclose(scaled.overall_contribution_change_from_v11.sum(),delta,rtol=0,atol=1e-12)
 np.testing.assert_allclose(np.mean([v for k,v in REPORTED.items() if k!='overall']),REPORTED['overall'],rtol=0,atol=1e-12)
 assert sum(REPORTED[k]>baseline['metrics'][k] for k in REPORTED if k!='overall')==4
 raw=pd.DataFrame([{'metric':k,'v11':baseline['raw_metrics'][k],'v13':api[k],'change_from_v11':api[k]-baseline['raw_metrics'][k],
  'relative_change':api[k]/baseline['raw_metrics'][k]-1,'improved_from_v11':bool(api[k]>baseline['raw_metrics'][k] if higher else api[k]<baseline['raw_metrics'][k])} for k,higher in RAW.items()])
 ready=json.loads((A/'submission_ready.json').read_text(encoding='utf-8'));audit=json.loads((A/'independent_validation.json').read_text(encoding='utf-8'))
 assert audit['passed'] and ready['sha256']==audit['package_sha256']==seal['package_sha256']
 assert (A.parent/'prediction.vcc').stat().st_size==ready['bytes'] and model['candidate']=='unknown_kd'
 for path in [ROOT/'analysis/leaderboard_result.json',*ROOT.glob('experiments/v*/analysis/leaderboard_result.json')]:
  if path.parent!=A:assert REPORTED['overall']>json.loads(path.read_text(encoding='utf-8'))['overall']
 result={'source':'User-reported CLI publication confirmed by read-only vcc --json status.','recorded_utc':datetime.now(timezone.utc).isoformat(),
  'entry_id':ENTRY,'model_name':api['model_name'],'final_status':'published','submitted_utc':api['submission_date'],
  'rank_at_user_report':558,'rank_at_api_check':api['rank'],'overall':REPORTED['overall'],'metrics':{k:v for k,v in REPORTED.items() if k!='overall'},
  'raw_metrics':{k:api[k] for k in RAW},'partition':api['partition'],'panel':api['panel_id'],'anchor_set':api['anchor_version'],
  'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'upload_duration_seconds':987,
  'overall_change_from_v11':delta,'relative_overall_change_from_v11':delta/baseline['overall'],
  'raw_relative_change':{k:api[k]/baseline['raw_metrics'][k]-1 for k in RAW},'release_manifest_sha256':seal_hash,
  'recorder_sha256':build.sha256(Path(__file__)),
  'interpretation':'Best recorded project validation score; modest improvement, four scaled components up, reach and overlap down. Local confirmation remains inconclusive. A single submission is not proof of statistical significance or final-round generalization.'}
 build.write_json(A/'published_api_status.json',api);build.write_json(A/'leaderboard_result.json',result)
 scaled.to_csv(A/'leaderboard_comparison.csv',index=False);raw.to_csv(A/'raw_metric_comparison.csv',index=False)
 archive=A/'generation_submission_document.md';expected=seal['files']['V13_SUBMISSION.md']
 if not archive.exists():
  assert build.sha256(ROOT/'V13_SUBMISSION.md')==expected;archive.write_bytes((ROOT/'V13_SUBMISSION.md').read_bytes())
 assert build.sha256(archive)==expected
 block=f"Published as `{ENTRY}`, overall **{REPORTED['overall']:.9f}**, a gain of **{delta:+.9f}** ({delta/baseline['overall']:.2%}) over V11 on the same panel and anchors. This is the best recorded project validation score. Rank was **558** in the user report and **{api['rank']}** at the API check; ranks change with other entries."
 original=archive.read_text(encoding='utf-8');title,rest=original.split('\n',1)
 rest=rest.replace('**Exploratory candidate: local tests did not confirm improvement over V11. V11 remains the recommended best-supported submission.**',
  '**Generation-time assessment:** local count confirmation was inconclusive because some primary NMAE panels were ineligible. The model was packaged as exploratory. This remains the historical validation outcome despite the later small leaderboard gain.')
 rest=rest.replace('It has not been uploaded.', 'It is now published; the original generation-time report is archived in `experiments/v13_response_programs/analysis/generation_submission_document.md`.')
 rest=rest.replace('- Published reference: V11, overall 0.137003974; no V13 leaderboard result exists.', '- Prior published reference: V11, overall 0.137003974. The V13 publication and metric comparison are recorded in [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).')
 rest=rest.replace('## Upload from Windows Command Prompt','## Historical upload command (already published)')
 rest=rest.replace("An upload is a leaderboard experiment. If interrupted after scoring starts, use `vcc status <entry-id> --wait` instead of creating another entry. V11's package and its published result remain intact.",f'Inspect the existing entry with `vcc status {ENTRY}`. The historical submit command creates a duplicate entry. V11 remains preserved for comparison.')
 (ROOT/'V13_SUBMISSION.md').write_text(title+'\n\n'+block+' See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).\n\n'+rest.lstrip('\n'),encoding='utf-8')
 report=f"""# V13 published leaderboard result

{block}

## Scaled comparison

{table(scaled)}
All displayed components are scaled scores: higher is better, including MSE and NMAE. Four components improve. Expression accuracy contributes +0.001322 to overall, partly offset by reach (-0.000741). Contributions are component changes divided by six; they reconcile exactly with the overall gain and do not establish causal attribution.

## Raw comparison

{table(raw)}
Raw expression error decreases **{-result['raw_relative_change']['expr_mse_unbiased_capped_norm']:.2%}**; raw reach decreases **{-result['raw_relative_change']['de_wilcoxon_direction_reach_raw']:.2%}**. Fidelity is still weak and its scaled component remains negative. These results support a small calibration improvement, not a broad solution to context transfer.

## Interpretation and next tests

V13 corrects unquantifiable on-target RNA being amplified as weak measured knockdown. It also includes small pooling-precision differences recorded before publication. Its population emitter and seed match V11. A single published result does not isolate these effects, establish statistical significance, or demonstrate performance on unseen final-round contexts. The frozen local count-confirmation result remains inconclusive; do not rewrite it as a pass.

Preserve V13 as the best recorded validation-score reference and V11 as a comparison with better reach and overlap. Prioritize a separately frozen test of the published construct-level K562 knockdown metadata, and coverage guards that preserve V11's reach. Neutral-population false-DE bias and unrealistic dispersion remain unresolved. Use held-family validation before another upload; this score alone does not justify tuning a larger model.

## Provenance and reproduction

Submitted `{api['submission_date']}`. Upload took 16 minutes 27 seconds in the user report. Partition `{api['partition']}`, panel `{api['panel_id']}`, anchors `{api['anchor_version']}`. Package SHA-256: `{ready['sha256']}`. The status API does not expose the uploaded hash; the user's upload record links the local audited artifact to this entry.

All frozen code, selection, confirmation, generation and integrity records remain unchanged, including generation-time `uploaded=false`. The original report is archived byte-for-byte. The release manifest still describes the generation-time files; `publication_manifest.json` maps its original `V13_SUBMISSION.md` hash to the archive and separately hashes current publication documents. Rerun this recorder after the original document generator to restore the publication overlay.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
.\\.competition_env\\Scripts\\python.exe record_v13_result.py --api-file experiments/v13_response_programs/analysis/published_api_status.json
```
"""
 (ROOT/'V13_LEADERBOARD_RESULTS.md').write_text(report,encoding='utf-8')
 readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8')
 text=text.replace('Nine submissions (first, V2, V4, V5, V6, V7, V9, V10 and V11) have been published.', 'Ten submissions (first, V2, V4, V5, V6, V7, V9, V10, V11 and V13) have been published.')
 text=text.replace('V11 is the best recorded validation submission, with overall score **0.137003974**, improving on V10 by **0.020249853**.', 'V13 is the best recorded validation submission, with overall score **0.137883134**, improving on V11 by **0.000879160** (0.64%).')
 text=text.replace('See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md) for the current result and [V11_SUBMISSION.md](V11_SUBMISSION.md) for reproducibility.', 'See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md) for the current result and [V13_SUBMISSION.md](V13_SUBMISSION.md) for reproducibility.')
 lines=text.splitlines()
 for i,line in enumerate(lines):
  if line.startswith('The latest and best published submission is'):lines[i]=f'The latest and best published submission is V13 (`{ENTRY}`), overall **{REPORTED["overall"]:.9f}**. See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).'
 text='\n'.join(lines)+'\n';text=text.replace('V11 is now the best published reference.', 'V13 is now the best published reference.').replace('V11 is now best.', 'V13 is now best.')
 marker='<!-- V13_RESPONSE_TRANSFER -->';assert marker in text;text=text[:text.index(marker)].rstrip()+'\n\n'+marker+'\n\n## V13 published result\n\n'+block+' Four components improve; reach and overlap decline. Local confirmation remains inconclusive. Official packaging and independent full-payload/provenance/replay checks passed. See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md), [V13_STRATEGY.md](V13_STRATEGY.md) and [V13_SUBMISSION.md](V13_SUBMISSION.md).\n'
 readme.write_text(text,encoding='utf-8')
 verify_generation(seal);assert build.sha256(A/'release_manifest.json')==seal_hash
 files=['published_api_status.json','leaderboard_result.json','leaderboard_comparison.csv','raw_metric_comparison.csv','generation_submission_document.md']
 build.write_json(A/'publication_manifest.json',{'entry_id':ENTRY,'recorder_sha256':build.sha256(Path(__file__)),
  'publication_artifact_sha256':{n:build.sha256(A/n) for n in files},'publication_document_sha256':{n:build.sha256(ROOT/n) for n in ['V13_SUBMISSION.md','V13_LEADERBOARD_RESULTS.md']},
  'release_manifest_sha256':seal_hash,'release_document_overlay':{'V13_SUBMISSION.md':'experiments/v13_response_programs/analysis/generation_submission_document.md'},
  'comparison_baseline_sha256':build.sha256(baseline_path),'all_generation_file_hashes_verified_before_and_after':True,
  'recorder_dependency_sha256':{n:build.sha256(ROOT/n) for n in ['build_submission.py']},
  'scope':'Publication overlay only; no change to frozen predictor, local acceptance, populations or package.'})
 print(scaled.to_string(index=False));print('V13 publication recorded; frozen generation evidence unchanged.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--api-file',type=Path);main(p.parse_args().api_file)
