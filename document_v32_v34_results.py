"""Document completed V32-V34 experiments without confusing raw and host scores."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import v32_methods as m
ROOT=m.ROOT
D32=m.A;D33=ROOT/'experiments/v33_tabular_response/analysis';D34=ROOT/'experiments/v34_residual_flow/analysis'

def read(p):return json.loads(p.read_text())
def count_assessment():
 p=read(D34/'plan.json');t=read(D34/'treatment_results.json');null=read(D34/'null_results.json');means=read(D34/'mean_results.json');baseline={(r['family'],r['arm']):r for r in t['panels']};changes={};checks={};g=p['raw_guards']
 fields=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
 for arm in ['flow_mean25','flow_shape25']:
  diff={};study={}
  for field in fields:
   values={}
   for r in t['panels']:
    if r['arm']!=arm:continue
    b=baseline[(r['family'],'v24')];a=r.get(field);bb=b.get(field)
    if a is None or bb is None:continue
    val=a/bb-1 if field=='common_expression_ratio' and bb>0 else a-bb
    values.setdefault(r['study'],[]).append(val)
   study[field]={k:float(np.mean(v)) for k,v in values.items()};diff[field]=float(np.mean(list(study[field].values()))) if len(study[field])==3 else None
  flags={'all_metrics_resolved':all(diff[k] is not None for k in fields)}
  limits={'pds':('pds_drop',1),'fidelity':('fidelity_drop',1),'reach':('reach_drop',1),'jaccard':('jaccard_drop',1),'nmae':('nmae_worsening',-1),'common_expression_ratio':('relative_expression_worsening',-1)}
  for k,(guard,sign) in limits.items():flags[k+'_guard']=diff[k] is not None and sign*diff[k]>=-g[guard]
  flags['family_fidelity_guard']=len(study['fidelity'])==3 and min(study['fidelity'].values())>=-g['family_fidelity_drop']
  flags['family_reach_guard']=len(study['reach'])==3 and min(study['reach'].values())>=-g['family_reach_drop']
  flags['two_DE_metrics_improve']=sum(diff[k] is not None and diff[k]>0 for k in ['fidelity','reach','jaccard'])>=2
  for r in null['records']:
   if r['arm']!=arm:continue
   flags['null_'+r['family']]=r['null_DE_calls']<=max(20,r['real_null_DE_calls']+10)
  flags['mean_stage_passed']=means['passes'];changes[arm]={'equal_study_changes':diff,'per_study_changes':study};checks[arm]={'checks':flags,'passes':all(flags.values())}
 result={'plan_sha256':m.sha(D34/'plan.json'),'changes':changes,'decisions':checks,'scope':'Fixed diagnostic subset, previously spent public biology. Raw score differences only; no host-anchor scaling. Shape arm descriptive and cannot override failed prespecified mean stage.'};m.write(D34/'count_assessment.json',result);return result

def main():
 for folder in [D32,D33,D34]:assert read(folder/'independent_audit.json')['passed']
 for path in [D33/'prediction_replay_audit.json',D34/'numpy_inference_audit.json',D34/'count_contract_audit.json']:assert read(path)['passed']
 r32=read(D32/'corrected_results.json');r33=read(D33/'results.json');r34=read(D34/'mean_results.json');count=count_assessment();rows=[('V32 paired response map',r32['summary']['paired'],r32['passes']),('V33 TabICLv2',r33['summary']['tabicl'],r33['passes']),('V34 conditional flow mean25',r34['summary']['flow_mean25'],r34['passes'])]
 lines=['# Three new model tests: V32–V34','', 'V24 remains the published comparator (overall 0.14299785521745675). All three hypotheses were implemented and tested in sequence. These are public development results, not leaderboard forecasts. The latest donor reserve was not opened.','', '| Primary model | Squared LFC error change vs V24 | Absolute LFC error change vs V24 | Mean gate |','|---|---:|---:|---|']
 for name,s,passed in rows:lines.append(f"| {name} | {(s['mse_ratio']-1)*100:+.3f}% | {(s['mae_ratio']-1)*100:+.3f}% | {'PASS' if passed else 'FAIL'} |")
 lines+=['','Negative error changes are improvements. V32/V33 use masked public summary-coordinate LFC errors; V34 uses approximate normalized mean-composition errors against real cells. Their percentages have different evaluation denominators and must not be ranked as interchangeable scores. Each compares against its matched held-family V24 baseline.','', '## V32: paired source-response transfer','',f"Corrected paired MSE ratio {r32['summary']['paired']['mse_ratio']:.9f}; identity-shuffled {r32['summary']['identity_shuffled']['mse_ratio']:.9f}; direct-copy {r32['summary']['source_copy']['mse_ratio']:.9f}.", '', 'The first implementation reapplied expression gating to V24 fallback. The separate correction preserves unsupported coordinates exactly and replays saved models with the same nested selection grid. Original results are retained for traceability. The corrected result, rather than the initial run, determines the conclusion.', '', 'All twelve ridge fits passed independently recomputed masked-gradient checks. All nine outer comparisons were replayed with independent float64 error accumulation; fitting/validation target separation and whole-family exclusions were checked.','', '## V33: tabular foundation model','', '| Matched algorithm, fixed alpha0.5 | Squared LFC error change | Absolute LFC error change |','|---|---:|---:|']
 for arm,s in r33['summary'].items():lines.append(f"| {arm} | {(s['mse_ratio']-1)*100:+.3f}% | {(s['mae_ratio']-1)*100:+.3f}% |")
 lines+=['','The pilot uses eight response programs, 512 route-balanced public in-context examples and two TabICL ensemble members. Ridge/CatBoost use identical inputs/labels. Their results are comparisons, not post-hoc alternative release selections. Checkpoint bytes, family/target boundaries, finite outputs and score aggregation were audited; every GPU forward was not independently reimplemented.','', '## V34: conditional flow','',f"Primary learned-mean CPM error ratio: {r34['summary']['flow_mean25']['cpm_mse_ratio']:.9f}.", '', 'Three independent outer-family flow models were trained for eight epochs each. Epoch selection used generated held-target source means, not outer scores or observed treated means supplied to the decoder. The architecture predicts 1,024 measured genes and retains V24 elsewhere. It does not receive the measured source-response tokens available to V24, so this comparison combines architecture and intervention-information differences. Missing STRING embeddings use a zero descriptor. Mean-only and flow-shape count diagnostics were run regardless of earlier model failures. The shape arm preserves the mean25 reference gene totals exactly, cell libraries exactly, and CPM mean L1 <=0.0005. After the initial integer projection failed in D1 (0.0007056), a uniform16-depth-bin allocation repaired the numerical contract; it also constrains learned gene-depth associations. Initial partial null results and failure evidence are retained.', '', '| Count diagnostic | PDS change | FID change | REACH change | JAC change | NMAE change | Relative expression error change |','|---|---:|---:|---:|---:|---:|---:|']
 for arm,rec in count['changes'].items():
  vals=rec['equal_study_changes'];lines.append('| '+arm+' | '+' | '.join('unresolved' if vals[k] is None else f'{vals[k]:+.6f}' for k in ['pds','fidelity','reach','jaccard','nmae','common_expression_ratio'])+' |')
 lines+=['','### Control-null calibration','', '| Context | Real NTC comparison DE calls | V24 emitter DE calls | Flow shape DE calls |','|---|---:|---:|---:|']
 null_rows=read(D34/'null_results.json')['records']
 for family in sorted({r['family'] for r in null_rows}):
  b=next(r for r in null_rows if r['family']==family and r['arm']=='v24');q=next(r for r in null_rows if r['family']==family and r['arm']=='flow_shape25');lines.append(f"| {family} | {b['real_null_DE_calls']} | {b['null_DE_calls']} | {q['null_DE_calls']} |")
 lines+=['','The actual NTC-vs-NTC comparisons have0–1 DE calls. Both synthetic emitters retain hundreds to thousands: the flow reduces this artifact only slightly. Null FID is not a treatment-recovery result. V29 previously showed that eliminating null calls alone also does not guarantee better treatment prediction.','', 'PDS/FID/REACH/JAC increases are favorable; NMAE/expression-error decreases are favorable. Undefined quantities remain null and fail resolution guards. The count subset consists of the first four cases in each of five contexts; D1 cases have 200 cells and are explicitly lower-power diagnostics. It does not estimate performance on the complete 300-target competition panel.','', '## Limits and decision','', 'This tests three specific, reproducible configurations. A negative pilot does not rule out its entire model class. There are few independent training lineages, most target/context examples come from reused public sources, and matched source/query assays differ. In particular, source-observed response transfer still must extrapolate to a recipient with only controls; recipient-anchored paper results do not remove that difficulty.','']
 if not any(v[2] for v in rows):lines.append('No primary mean candidate passed. No new submission package was generated, and V24 remains the preferred upload. Further tuning on these same outer outcomes would be exploratory and would require new confirmation evidence.')
 else:lines.append('A primary mean gate passed; count/null outcomes and the recorded scope still govern whether it can advance. No automatic upload occurred.')
 lines+=['','See [V32_V34_STRATEGY.md](V32_V34_STRATEGY.md) for inputs, exact commands, corrections and literature. Plans, score records, checkpoints receipts and independent audit results are retained under their respective experiment analysis directories.','']
 lines+=['![Mean-error comparison](experiments/v34_residual_flow/analysis/comparison.png)','']
 (ROOT/'V32_V34_RESULTS.md').write_text('\n'.join(lines),encoding='utf-8');fig,ax=plt.subplots(figsize=(9,4.2));names=[v[0] for v in rows];x=np.arange(3);ax.bar(x-.18,[(v[1]['mse_ratio']-1)*100 for v in rows],.36,label='Squared LFC error');ax.bar(x+.18,[(v[1]['mae_ratio']-1)*100 for v in rows],.36,label='Absolute LFC error');ax.axhline(0,color='black',lw=.8);ax.set_xticks(x,names,fontsize=9);ax.set_ylabel('Error change versus matched V24 (%)\nLower is better');ax.legend();ax.set_title('Three prespecified public-development pilots');fig.text(.5,.015,'V32/V33: summary-coordinate errors; V34: real-cell mean composition. Different denominators.',ha='center',fontsize=8);fig.tight_layout(rect=(0,.04,1,1));fig.savefig(D34/'comparison.png',dpi=180);plt.close(fig);m.write(D34/'three_model_decision.json',{'V32_passes':r32['passes'],'V33_passes':r33['passes'],'V34_mean_passes':r34['passes'],'V34_count_decisions':count['decisions'],'new_submission_generated':False,'reserve_opened':False,'published_champion':'V24','published_overall':.14299785521745675});print('\n'.join(lines[:14]))
if __name__=='__main__':main()