"""Freeze the pre-outcome V35/V36 release comparison; never imply private gains."""
import json
import numpy as np
import v35_methods as m

RULE={'mean_selection':'Choose v35 or v36 with lowest equal-study squared LFC error on the44 frozen raw cases, among candidates passing count noninferiority guards if any. If none pass, package the lowest-mean-error candidate as exploratory only. Never change the V24 champion record.','dispersion_selection':'For a selected v36 mean, adopt its dispersion arm only if it passes all count guards and improves at least2 of FID,REACH,JAC,NMAE with no expression-error regression versus the unchanged-emitter v36 arm. Otherwise unchanged V24 emitter.','guards':{'maximum_relative_expression_error':1.02,'maximum_FID_drop':.02,'maximum_PDS_drop':.02,'null_DE_calls':'no more than matched V24 for every context'},'reserve':'Never open D3 for exploratory packaging; no blind improvement claim.','submission':'Prepare validated upload package; no automatic upload, existing packages unchanged.'}

def freeze():
 p=m.P36/'release_comparison_plan.json';assert not p.exists();m.write(p,{'parent_plan_sha256':m.sha(m.A/'plan.json'),'code_sha256':m.sha(m.ROOT/'select_v35_v36_release.py'),'rule':RULE});print('Release comparison frozen',flush=True)

def select():
 plan=json.loads((m.P36/'release_comparison_plan.json').read_text());assert plan['code_sha256']==m.sha(m.ROOT/'select_v35_v36_release.py')
 means=json.loads((m.P36/'mean_results.json').read_text());counts=json.loads((m.P36/'treatment_results.json').read_text());null=json.loads((m.P36/'null_results.json').read_text());summaries={}
 metrics=['fidelity','reach','jaccard','nmae','pds','common_expression_ratio']
 for arm in ['v24','v35','v36','v36_dispersion']:
  studies={}
  for study in m.OUTERS:
   rr=[r for r in counts['panels'] if r['study']==study and r['arm']==arm];studies[study]={}
   for k in metrics:
    vals=[r[k] for r in rr if isinstance(r.get(k),(float,int)) and np.isfinite(r[k])];studies[study][k]=float(np.mean(vals)) if vals else None
  summaries[arm]={k:float(np.mean([v[k] for v in studies.values()])) if all(v[k] is not None for v in studies.values()) else None for k in metrics}
 def guard(arm):
  s=summaries[arm];b=summaries['v24'];checks={'expression':s['common_expression_ratio'] is not None and s['common_expression_ratio']<=1.02*b['common_expression_ratio'],'FID':s['fidelity'] is not None and s['fidelity']>=b['fidelity']-.02,'PDS':s['pds'] is not None and s['pds']>=b['pds']-.02,'null':all(r['null_DE_calls']<=next(v['null_DE_calls'] for v in null['records'] if v['family']==r['family'] and v['arm']=='v24') for r in null['records'] if r['arm']==arm)}
  return {'checks':checks,'passes':all(checks.values())}
 guards={a:guard(a) for a in ['v35','v36','v36_dispersion']};eligible=[a for a in ['v35','v36'] if guards[a]['passes']];chosen=min(eligible or ['v35','v36'],key=lambda a:means['summary'][a]['mse_ratio']);arm=chosen
 if chosen=='v36':
  s=summaries['v36_dispersion'];b=summaries['v36'];better=sum(s[k] is not None and b[k] is not None and (s[k]<b[k] if k=='nmae' else s[k]>b[k]) for k in ['fidelity','reach','jaccard','nmae'])
  if guards['v36_dispersion']['passes'] and better>=2 and s['common_expression_ratio']<=b['common_expression_ratio']:arm='v36_dispersion'
 summary=json.loads(((m.A if chosen=='v35' else m.P36)/'summary_results.json').read_text());strict=summary['passes'] and means['summary'][chosen]['mse_ratio']<=.99 and means['summary'][chosen]['mae_ratio']<=.995 and guards[arm]['passes']
 out={'selected_arm':arm,'model_version':35 if chosen=='v35' else 36,'strict_public_guards_passed':bool(strict),'status':'public-development candidate; unconfirmed on new biology' if strict else 'exploratory candidate; strict improvement guards not all passed','confirmed_improvement':False,'available_public_guards_passed':bool(strict),'count_summaries':summaries,'count_guards':guards,'raw_mean_summaries':means['summary'],'selection_plan_sha256':m.sha(m.P36/'release_comparison_plan.json'),'V24_champion_unchanged':True,'reserve_opened':False,'uploaded':False,'scope':'Raw public metrics are not private scaled VCC scores. Package requested by user; absence of strict gains is reported, not converted into a promise.'}
 m.write(m.P36/'selection.json',out);print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','select']);args=p.parse_args();freeze() if args.stage=='freeze' else select()