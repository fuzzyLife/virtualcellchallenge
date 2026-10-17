"""Select one new version before additional count-confirmation outcomes are opened."""
import argparse
import build_submission as b
import v46_population_mixture as one
import v47_kolf_consensus as two
A=two.A
def freeze():
 assert not (A/'cross_version_rule.json').exists()
 for m in [one,two]:assert not (m.A/'confirmation_results.json').exists()
 b.write_json(A/'cross_version_rule.json',{'selector_sha256':b.sha256(b.ROOT/'select_v46_v47_final.py'),'plan_sha256':{'46':b.sha256(one.A/'plan.json'),'47':b.sha256(two.A/'plan.json')},'rule':'Only candidates passing every own frozen development/mean/control guard qualify. Highest oriented six-metric relative gain wins; if within.001 prefer V47 independent experimental corroboration. Lock only one version before opening any remaining24 count case. Do not reselect on confirmation outcomes.','known_outcomes':'V38–V45 complete development and partial V46 development visible; V47 predictions not yet completed. Public biology repeatedly used. Additional24 count outcomes and D3 unopened.','private_score_known':False})
 print('Cross-version selection rule frozen',flush=True)
def select():
 rule=two.read(A/'cross_version_rule.json');assert b.sha256(b.ROOT/'select_v46_v47_final.py')==rule['selector_sha256'];candidates=[]
 for v,m in [(46,one),(47,two)]:
  assert not (m.A/'confirmation_results.json').exists();assert b.sha256(m.A/'plan.json')==rule['plan_sha256'][str(v)]
  s=m.read(m.A/'selection.json')
  for arm in s['eligible']:candidates.append({'version':v,'arm':arm,'gain':s['relative_gain'][arm],'selection_sha256':b.sha256(m.A/'selection.json')})
 if candidates:
  best=max(candidates,key=lambda x:x['gain']);near=[c for c in candidates if best['gain']-c['gain']<.001]
  if any(c['version']==47 for c in near):best=max([c for c in near if c['version']==47],key=lambda x:x['gain'])
 else:best=None
 b.write_json(A/'locked_cross_version.json',{'selected':best,'candidates':candidates,'rule_sha256':b.sha256(A/'cross_version_rule.json'),'additional_count_confirmation_opened':False,'D3_opened':False})
 print('Locked cross-version finalist',best,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','select']);a=p.parse_args();globals()[a.stage]()

