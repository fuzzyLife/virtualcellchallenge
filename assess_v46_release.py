"""Lock a passing public-development arm, then assess untouched count confirmation."""
import argparse
import numpy as np
import v46_population_mixture as m
from run_v46_count_grid import METRICS
from run_v37_ensemble import aggregate
def guards(s,base):
 return {'resolved':all(s.get(k) is not None and base.get(k) is not None for k in METRICS),
 'pds':s['pds']>=base['pds']-.01,'expression':s['common_expression_ratio']<=base['common_expression_ratio']*1.02,
 'fid':s['fidelity']>=base['fidelity']-.01,'reach':s['reach']>=base['reach']-.005,
 'jac':s['jaccard']>=base['jaccard']-.001,'nmae':s['nmae']<=base['nmae']*1.01}
def score(s,base):
 return float(np.mean([(s[k]-base[k])/max(abs(base[k]),.01)*(-1 if k in ('nmae','common_expression_ratio') else 1) for k in METRICS]))
def lock():
 m.check();s=m.read(m.A/'selection.json');assert s['selected']!='v24' and s['selected'] in s['eligible']
 assert not (m.A/'locked_finalist.json').exists()
 m.write(m.A/'locked_finalist.json',{'selected':s['selected'],'selection_sha256':m.b.sha256(m.A/'selection.json'),'plan_sha256':m.b.sha256(m.A/'plan.json'),'confirmation_was_unopened':not (m.A/'confirmation_results.json').exists(),'D3_opened':False,'scope':'Candidate fixed before remaining24 count cases at2 seeds; biological means already used publicly.'})
 print('V46 locked',s['selected'],flush=True)
def assess():
 m.check();lock=m.read(m.A/'locked_finalist.json');assert m.b.sha256(m.A/'selection.json')==lock['selection_sha256']
 arm=lock['selected'];summary=m.read(m.A/'confirmation_summary.json');base=summary['v24'];s=summary[arm];tests=guards(s,base);tests['positive_average_gain']=score(s,base)>0
 tests['two_metrics_improve']=sum(s[k]<base[k] if k in ('nmae','common_expression_ratio') else s[k]>base[k] for k in METRICS)>=2
 panels=m.read(m.A/'confirmation_results.json')['panels'];study_checks={}
 for study in ['H1','Orion','CD4']:
  vv=[p for p in panels if p['study']==study];vals={a:{k:float(np.mean([np.mean([p[k] for p in vv if p['arm']==a and p['family']==f]) for f in sorted({p['family'] for p in vv})])) for k in ['pds','fidelity','reach']} for a in ['v24',arm]}
  study_checks[study]={k:vals[arm][k]>=vals['v24'][k]-.03 for k in vals[arm]}
 tests['worst_study']=all(all(x.values()) for x in study_checks.values())
 passed=all(tests.values())
 m.write(m.A/'release_decision.json',{'release':passed,'selected':arm,'all_public_guards_passed':passed,'confirmation_checks':tests,'study_checks':study_checks,'oriented_relative_gain':score(s,base),'confirmation_summary':summary,'locked_finalist_sha256':m.b.sha256(m.A/'locked_finalist.json'),'confirmation_results_sha256':m.b.sha256(m.A/'confirmation_results.json'),'private_score_known':False,'D3_opened':False,'assessor_sha256':m.b.sha256(m.ROOT/'assess_v46_release.py')})
 print('V46 public confirmation',passed,tests,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['lock','assess']);a=p.parse_args();globals()[a.stage]()




