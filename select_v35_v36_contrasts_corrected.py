"""Repair selector handling of undefined per-context metrics; same thresholds.

HCT116 has too few reference DE genes for NMAE in every arm. Match the original
primary selector: exclude undefined contexts within a study; if the whole study
is unresolved, mark the aggregate unresolved and fail its guard. Never use zero.
"""
import json,numpy as np
import v35_methods as m
from run_v35_v36_contrasts import A,ARMS

def aggregate(rows,arm,metric):
 studies=[]
 for study in m.OUTERS:
  vals=[r[metric] for r in rows if r['arm']==arm and r['study']==study and isinstance(r.get(metric),(float,int)) and np.isfinite(r[metric])]
  if not vals:return None
  studies.append(float(np.mean(vals)))
 return float(np.mean(studies))

def main():
 original=m.P36/'precontrast_selection.json';selection=json.loads((m.P36/'selection.json').read_text());assert 'contrastive_correction' not in selection
 if original.exists():assert original.read_bytes()==(m.P36/'selection.json').read_bytes()
 else:original.write_bytes((m.P36/'selection.json').read_bytes())
 result=json.loads((A/'results.json').read_text());base=selection['count_summaries']['v24'];metrics=['fidelity','reach','jaccard','nmae','pds','common_expression_ratio'];summaries={};checks={};eligible=[]
 for arm in ARMS:
  s={k:aggregate(result['panels'],arm,k) for k in metrics};summaries[arm]=s;mean=result['mean_summary'][arm];resolved=all(s[k] is not None for k in metrics);guard=resolved and s['common_expression_ratio']<=1.02*base['common_expression_ratio'] and s['fidelity']>=base['fidelity']-.02 and s['pds']>=base['pds']-.02 and s['nmae']<=1.01*base['nmae'];improved=mean['mse_ratio']<1 or mean['mae_ratio']<1 or any(s[k] is not None and s[k]>base[k] for k in ['fidelity','reach','jaccard']);checks[arm]={'metrics_resolved':resolved,'count_guard':bool(guard),'mean_noninferiority':mean['mse_ratio']<=1.002 and mean['mae_ratio']<=1.002,'some_gain':bool(improved)}
  if all(checks[arm].values()):eligible.append(arm)
 current=selection['selected_arm'];c=selection['count_summaries'][current];cm=selection['raw_mean_summaries']['v35' if selection['model_version']==35 else 'v36'];current_ok=c['common_expression_ratio']<=1.02*base['common_expression_ratio'] and c['fidelity']>=base['fidelity']-.02 and c['pds']>=base['pds']-.02 and c['nmae']<=1.01*base['nmae'] and cm['mse_ratio']<=1.002 and cm['mae_ratio']<=1.002
 if current_ok:eligible.append(current)
 choice=min(eligible,key=lambda a:result['mean_summary'][a]['mse_ratio'] if a in ARMS else cm['mse_ratio']) if eligible else current
 selection['count_summaries'].update(summaries);selection['contrastive_correction']={'plan_sha256':m.sha(A/'plan.json'),'results_sha256':m.sha(A/'results.json'),'prior_selection_sha256':m.sha(original),'checks':checks,'eligible':eligible,'selected':choice,'exploratory_after_primary_outcomes':True,'selector_correction_sha256':m.sha(m.ROOT/'select_v35_v36_contrasts_corrected.py'),'undefined_context_rule':'Same primary aggregation: omit per-context unresolved NMAE, require every biological study resolved. HCT116 NMAE undefined in all arms, never replaced by zero.'};selection['raw_mean_summaries'].update(result['mean_summary'])
 if choice in ARMS:selection['selected_arm']=choice;selection['model_version']=35;selection.pop('dispersion_model_version',None)
 m.write(m.P36/'selection.json',selection);m.write(A/'decision.json',selection['contrastive_correction']);print(json.dumps({'selected':choice,'checks':checks,'means':result['mean_summary'],'counts':summaries},indent=2),flush=True)
if __name__=='__main__':main()