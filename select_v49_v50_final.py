"""Select one V49/V50 candidate before opening additional count confirmation."""
import argparse
import build_submission as b
import v49_balanced_population as one
import v50_response_prior as two
A=two.A
def freeze():
 assert not (A/'cross_version_rule.json').exists()
 for m in [one,two]:assert not (m.A/'confirmation_results.json').exists()
 b.write_json(A/'cross_version_rule.json',{'selector_sha256':b.sha256(b.ROOT/'select_v49_v50_final.py'),'plan_sha256':{'49':b.sha256(one.A/'plan.json'),'50':b.sha256(two.A/'plan.json')},'rule':'Only own fully passing development candidates qualify. Highest oriented six-raw-metric relative gain wins; within.001 prefer V50 mobility_prior if qualified, otherwise highest gain. Lock only one version/arm before additional24 count cases at2 seeds. No reselection after confirmation.','known_outcomes':'V38–V48 known, V49 running; V50 randomized checkpoint development already known from V47 negative control, averaged prior not yet count-tested. Explicit outcome-informed new interpretation. Additional24 count cases and D3 unopened.','private_score_known':False})
 print('V49/V50 cross-version rule frozen',flush=True)
def select():
 rule=two.read(A/'cross_version_rule.json');assert b.sha256(b.ROOT/'select_v49_v50_final.py')==rule['selector_sha256'];candidates=[]
 for v,m in [(49,one),(50,two)]:
  assert not (m.A/'confirmation_results.json').exists();assert b.sha256(m.A/'plan.json')==rule['plan_sha256'][str(v)];s=m.read(m.A/'selection.json')
  for a in s['eligible']:candidates.append({'version':v,'arm':a,'gain':s['relative_gain'][a],'selection_sha256':b.sha256(m.A/'selection.json')})
 best=max(candidates,key=lambda x:x['gain']) if candidates else None
 if best:
  avg=[c for c in candidates if c['version']==50 and c['arm']=='mobility_prior' and best['gain']-c['gain']<.001]
  if avg:best=avg[0]
 b.write_json(A/'locked_cross_version.json',{'selected':best,'candidates':candidates,'rule_sha256':b.sha256(A/'cross_version_rule.json'),'additional_count_confirmation_opened':False,'D3_opened':False});print('Locked finalist',best,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','select']);a=p.parse_args();globals()[a.stage]()

