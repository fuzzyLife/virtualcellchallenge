"""Report study/seed stability from existing raw confirmation statistics."""
import numpy as np
import v50_response_prior as m
from assess_v50_release import METRICS,study_value,score
from run_v37_ensemble import aggregate
def main():
 result=m.read(m.A/'confirmation_results.json');panels=result['panels'];arm=m.read(m.A/'locked_finalist.json')['selected'];studies={}
 for study in ['H1','Orion','CD4']:
  rows=[p for p in panels if p['study']==study];studies[study]={a:{k:study_value(rows,a,k) for k in METRICS} for a in ['v24',arm]}
 seeds={}
 for offset in sorted({p['seed_offset'] for p in panels}):
  rows=[p for p in panels if p['seed_offset']==offset];s={a:{k:aggregate(rows,a,k) for k in METRICS} for a in ['v24',arm]};seeds[str(offset)]={'summary':s,'oriented_relative_gain':score(s[arm],s['v24'])}
 m.write(m.A/'confirmation_stability.json',{'studies':studies,'seeds':seeds,'scope':'Descriptive raw public diagnostics. Three reused study families and two generation seeds do not establish statistical significance or private leaderboard improvement. Same initial RNG seed does not guarantee identical Poisson draws when rates differ.','script_sha256':m.b.sha256(m.ROOT/'summarize_v50_confirmation.py')})
 print('V50 per-seed gains',{k:v['oriented_relative_gain'] for k,v in seeds.items()},flush=True)
if __name__=='__main__':main()

