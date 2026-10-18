"""Record official V55 publication separately from frozen pre-upload evidence."""
from pathlib import Path
import argparse,json,hashlib,math
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v55_mean50_diagnostic/analysis'
METRICS=('pds','mse','nmae','fid','reach','jac')
RAW=('pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard')
ENTRY='CTdmP6ishgalaMMz4XRe'

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scores(d):return {'overall':d['overall'],**d['metrics']}
def render():
 api_path=A/'published_api_status.json';api=read(api_path)
 assert api['entry_id']==ENTRY and api['status']=='published' and api['is_terminal']
 assert api['model_name']=='Fixed V24 V36 half-blend diagnostic v55'
 current={'overall':api['score_avg'],**{k:api['score_'+k] for k in METRICS}}
 assert all(math.isfinite(v) for v in current.values())
 assert abs(sum(current[k] for k in METRICS)/6-current['overall'])<1e-12
 bases={v:read(ROOT/f'experiments/{directory}/analysis/leaderboard_result.json') for v,directory in [('v24','v24_cd4_response'),('v36','v36_population_views')]}
 panel=(api['partition'],api['panel_id'],api['anchor_version'])
 for old in bases.values():assert panel==(old['partition'],old['panel'],old['anchor_set'])
 ready=read(A/'submission_ready.json');manifest=read(A/'final_manifest.json')
 assert ready['official_prep_passed'] and ready['sha256']==manifest['package_sha256']
 assert ready['sha256']=='b877efdd31565ed7a58bcf8f26384a8837e546221fab5733ab0602e5db713c1b'
 assert (A.parent/'prediction.vcc').stat().st_size==ready['bytes']==manifest['package_bytes']
 sealed={}
 for n,h in manifest['receipts_sha256'].items():assert sha(A/n)==h,n;sealed['analysis/'+n]=h
 for n,h in manifest['documentation_sha256'].items():assert sha(ROOT/n)==h,n;sealed[n]=h
 comparisons={v:{k:current[k]-scores(d)[k] for k in current} for v,d in bases.items()}
 for changes in comparisons.values():assert abs(sum(changes[k]/6 for k in METRICS)-changes['overall'])<1e-12
 contribution={k:comparisons['v24'][k]/6 for k in METRICS}
 reach_fraction=contribution['reach']/comparisons['v24']['overall']
 assessment=read(A/'assessment.json');public=assessment['summaries'];local_fid=public['mean50']['fidelity']-public['v24']['fidelity'];local_reach=public['mean50']['reach']-public['v24']['reach']
 assert local_fid>0 and comparisons['v24']['fid']<0
 assert local_reach<0 and comparisons['v24']['reach']<0
 result={'source':'Authenticated read-only official vcc --json status after exact model discovery on the public leaderboard.','review_date_client_timezone':'2026-10-18 Europe/Oslo','entry_id':ENTRY,'model_name':api['model_name'],'final_status':api['status'],'submitted_utc':api['submission_date'],'rank_at_api_check':api['rank'],'overall':current['overall'],'metrics':{k:current[k] for k in METRICS},'raw_metrics':{k:api[k] for k in RAW},'partition':panel[0],'panel':panel[1],'anchor_set':panel[2],'all_scaled_metrics_higher_is_better':True,'package_sha256':ready['sha256'],'upload_bytes':ready['bytes'],'package_link_provenance':'User confirmed completion after the supplied exact command; API independently confirms owned entry and model name but does not expose a remote upload checksum. Package hash is from preserved local validation receipts.','changes_from_v24':comparisons['v24'],'changes_from_v36':comparisons['v36'],'overall_contributions_vs_v24':contribution,'reach_fraction_of_net_overall_decline':reach_fraction,'best_recorded_project_validation_score':False,'preferred_submission':'experiments/v24_cd4_response/prediction.vcc','public_private_comparison':{'public_raw_fidelity_delta':local_fid,'private_scaled_fidelity_delta':comparisons['v24']['fid'],'public_raw_reach_delta':local_reach,'private_scaled_reach_delta':comparisons['v24']['reach'],'scope':'Signs can be compared; raw public and anchored private magnitudes have different units. Public panel/scoring differences remain.'},'api_receipt_sha256':sha(api_path),'frozen_final_manifest_sha256':sha(A/'final_manifest.json'),'sealed_preupload_sha256':sealed,'recorder_sha256':sha(Path(__file__)),'interpretation':'V55 loses overall to both V24 and V36. Only expression accuracy improves versus V24; REACH accounts for most of the net decline. Local fidelity gains fail to transfer, weakening the previous small-sample directional agreement. The historical strict V37 REACH rejection remains preserved. No significance or final-test generalization claim.'}
 user=read(A/'published_user_report.json') if (A/'published_user_report.json').exists() else None
 if user:
  assert user['model_name']==api['model_name'] and round(current['overall'],4)==user['rounded_overall']
  assert 1<=user['rank_at_user_report']<=user['teams_at_user_report']
  result.update(rank_at_user_report=user['rank_at_user_report'],teams_at_user_report=user['teams_at_user_report'],user_report_sha256=sha(A/'published_user_report.json'))
 lines=['# V55 published result','',
 f"Reviewed 18 October 2026 (Europe/Oslo). Official authenticated API confirms entry `{ENTRY}`, model `{api['model_name']}`, **published**, rank **{api['rank']} at this check**. Historical ranks change with other submissions; compare scores on the same panel and anchors.",'',
 f"Overall **{current['overall']:.9f}**, change versus V24 **{comparisons['v24']['overall']:+.9f}** ({comparisons['v24']['overall']/bases['v24']['overall']:+.2%}), versus V36 **{comparisons['v36']['overall']:+.9f}**. **V24 remains the best published project submission.**",'',
 '| Scaled metric, higher is better | V24 | V36 | V55 | V55 minus V24 | Overall contribution vs V24 |','|---|---:|---:|---:|---:|---:|']
 for k in current:
  lines.append(f"| {k} | {scores(bases['v24'])[k]:.9f} | {scores(bases['v36'])[k]:.9f} | {current[k]:.9f} | {comparisons['v24'][k]:+.9f} | "+('n/a' if k=='overall' else f'{contribution[k]:+.9f}')+' |')
 lines+=['',f"REACH contributes {contribution['reach']:+.9f} to overall, **{reach_fraction:.1%} of the net decline**. Expression accuracy improves, while the other five scaled components decline versus V24. This is arithmetic decomposition of the mean of six scores, not a causal attribution. MSE/NMAE columns above are scaled accuracies and increase when favorable; their raw errors decrease when favorable.",'',
 '## What the diagnostic establishes','',
 'The fixed half-blend failed to improve on either parent. Moving halfway between response effects did not interpolate the hosted score between the parent scores. Count normalization, rounding and DE ranking/significance are nonlinear; this aggregate outcome does not identify which operation or genes caused the decline.','',
 f"Public raw fidelity improved by {local_fid:+.9f}, while private scaled fidelity changed by {comparisons['v24']['fid']:+.9f}. The apparent4/4 local-to-official fidelity sign agreement from V54 now has a counterexample. It was a small, heterogeneous retrospective sample, not a reliable model-selection rule. Expression improved in both public and private comparisons, but this alone did not select a better overall model.",'',
 f"Public raw REACH changed by {local_reach:+.9f}; private scaled REACH changed by {comparisons['v24']['reach']:+.9f}. The direction of the loss was visible in the new public tests and in the original V37 rejection. Poor correlation of local REACH gains with private gains was not evidence that a REACH loss could safely be disregarded. Raw public/private scaled magnitudes cannot be equated.",'',
 'Retain V24 as reference. This diagnostic supplies no evidence for another similar residual blend. Before choosing another upload, test the mean-response versus distribution/rank bottleneck on public data with fixed matched comparisons and independent reserve discipline. No new training, package, reupload or reserve opening is part of this result check.','',
 '## Provenance and preserved evidence','',
 f"Partition `{panel[0]}`, panel `{panel[1]}`, anchors `{panel[2]}` match V24 and V36. API submission timestamp: `{api['submission_date']}` (UTC). Local package SHA-256 `{ready['sha256']}`; {ready['bytes']:,} bytes. The API does not expose an uploaded-file checksum; it confirms the owned model/entry, not remote byte identity.",'',
 'The original [V55 generation report](V55_MEAN50_RESULTS.md), protocol, ready-file receipt, final manifest, arrays and package remain unchanged, including their generation-time unknown-score/uploaded-false fields. This publication record is separate. The [official API receipt](experiments/v55_mean50_diagnostic/analysis/published_api_status.json) and [canonical result](experiments/v55_mean50_diagnostic/analysis/leaderboard_result.json) record the subsequent outcome.','',
 '~~~bat','cd /d L:\\promec\\Animesh\\Download\\Virtual',f'.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}',
 '"C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe" record_v55_result.py --check','~~~','']
 if user:lines.insert(4,f"The user-pasted page shows rank **{user['rank_at_user_report']}/{user['teams_at_user_report']}**; the authenticated API snapshot shows **{api['rank']}**. Both rank snapshots are retained. Score values match after rounding; no claim is made that these ranks are simultaneous.")
 return result,('\n'.join(lines)).encode('utf-8')

def main(check=False,amend_user_report=False):
 result,report=render();p=A/'leaderboard_result.json';doc=ROOT/'V55_LEADERBOARD_RESULTS.md'
 if check:
  assert read(p)==result and doc.read_bytes()==report
  print('PASS: published API scores, same anchors, mean-six arithmetic, both parent comparisons, all sealed receipts and publication report.')
  return
 if p.exists() and read(p)!=result and amend_user_report:
  old=read(p);allowed={'recorder_sha256','rank_at_user_report','teams_at_user_report','user_report_sha256'}
  assert {k:v for k,v in old.items() if k not in allowed}=={k:v for k,v in result.items() if k not in allowed},'Only recorder/user-rank provenance may change.'
  backup=A/'api_only_initial_leaderboard_result.json';assert not backup.exists();backup.write_bytes(p.read_bytes())
  olddoc=A/'api_only_initial_leaderboard_report.md';assert not olddoc.exists();olddoc.write_bytes(doc.read_bytes())
  p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n');doc.write_bytes(report)
 if p.exists():assert read(p)==result,'Refuse to overwrite a different published receipt.'
 else:p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 if doc.exists():assert doc.read_bytes()==report,'Refuse to overwrite a different publication report.'
 else:doc.write_bytes(report)
 f=ROOT/'README.md';body=f.read_text(encoding='utf-8-sig');start='<!-- v55-mean50:begin -->';end='<!-- v55-mean50:end -->';a=body.index(start);z=body.index(end,a)+len(end)
 block=start+'\n\n## V55 published diagnostic result\n\nV55 is published as `'+ENTRY+'`, overall **'+f"{result['overall']:.9f}"+'**, change **'+f"{result['changes_from_v24']['overall']:+.9f}"+'** versus V24. User page rank568/1568; API snapshot rank575 at the 18 October 2026 check. Both snapshots are retained. Only expression accuracy improves; REACH accounts for '+f"{result['reach_fraction_of_net_overall_decline']:.1%}"+' of the net overall decline. The public fidelity gain fails to transfer. **V24 remains the best published project submission.** [Published scores and interpretation](V55_LEADERBOARD_RESULTS.md), [preserved pre-upload report](V55_MEAN50_RESULTS.md), [protocol](V55_MEAN50_STRATEGY.md). All generation-time receipts and the package remain unchanged.\n\n'+end
 f.write_text(body[:a]+block+body[z:],encoding='utf-8',newline='\r\n')
 print(json.dumps({'overall':result['overall'],'rank':result['rank_at_api_check'],'delta_v24':result['changes_from_v24']['overall'],'delta_v36':result['changes_from_v36']['overall'],'reach_fraction_of_decline':result['reach_fraction_of_net_overall_decline']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--amend-user-report',action='store_true');args=p.parse_args();main(args.check,args.amend_user_report)
