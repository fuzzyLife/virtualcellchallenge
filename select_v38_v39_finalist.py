"""Lock one V38/V39 finalist before opening additional raw count cases."""
import json
from pathlib import Path
import v38_response_hurdle as m
ROOT=m.ROOT
def main():
 out=ROOT/'experiments/v38_response_hurdle/analysis'
 rule=out/'cross_version_rule.json'
 if not rule.exists():
  m.write(rule,{'script_sha256':m.b.sha256(ROOT/'select_v38_v39_finalist.py'),'known_outcomes':'All V38 means, source validation and zero-head ablation known; V38 development is still underway. V39 source models trained; V39 public mean and count results not yet inspected.','rule':'Require each version own frozen development/null guards and selected mean arm eligible. Among qualifying versions, retain V38 unless V39 average six-metric relative gain exceeds V38 by0.001; if only one qualifies choose it. Lock before additional24 raw count cases. Confirm original frozen count guards for that candidate across both seeds, plus each study PDS/FID/REACH no worse than V24 by0.03. No reselection after confirmation.'})
 assert m.read(rule)['script_sha256']==m.b.sha256(ROOT/'select_v38_v39_finalist.py')
 versions={}
 for v,folder in [(38,'v38_response_hurdle'),(39,'v39_kolf_support')]:
  a=ROOT/'experiments'/folder/'analysis';sel=m.read(a/'selection.json');mean=m.read(a/'mean_results.json');arm=sel['selected']
  if arm!='v24' and arm.removesuffix('_transport') in mean['eligible']:
   versions[v]={'version':v,'folder':folder,'selected':arm,'gain':sel['relative_gain'][arm],'selection_sha256':m.b.sha256(a/'selection.json')}
 selected=versions.get(38)
 if 39 in versions and (selected is None or versions[39]['gain']>selected['gain']+.001):selected=versions[39]
 p=out/'locked_finalist.json';assert not p.exists();m.write(p,{'selected':selected,'eligible':versions,'rule_sha256':m.b.sha256(rule),'confirmation_opened':False,'official_score_unknown':True})
 print(json.dumps(m.read(p),indent=2),flush=True)
if __name__=='__main__':main()

