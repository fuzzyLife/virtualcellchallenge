"""Apply V89's prespecified strong-support and bounded experimental-upload decision."""
import numpy as np
import v89_projection_release as m
from v57_pie import read,write,sha

def main():
 m.check();broad=read(m.A/'mean_results.json');summary=broad['summary']['interval25'];broad_safety=dict(profile=summary['profile_ratio']<=1.01,worst_study=max(x['profile_ratio'] for x in broad['studies']['interval25'].values())<=1.03,pds=summary['pds_change']>=-.01);seeds={};strong=broad['passes'];safe=all(broad_safety.values())
 for offset in [0,620000]:
  r=read(m.A/f'count_assessment_{offset}.json');base=r['summary']['v66'];v=r['summary']['interval25'];checks={}
  for k,slack,smaller in [('profile_MSE',.005,True),('pds',.005,False),('fidelity',.002,False),('jaccard',.0005,False),('nmae',.01,True)]:
   checks[k]=v[k] is not None and base[k] is not None and (v[k]<=(1+slack)*base[k] if smaller else v[k]>=base[k]-slack)
  improved=[k for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae'] if v[k] is not None and base[k] is not None and (v[k]<base[k] if k in ['profile_MSE','nmae'] else v[k]>base[k])];checks['two_gains']=len(improved)>=2;safe&=all(checks.values());strong&=r['passes'];seeds[str(offset)]=dict(summary=r['summary'],experimental_checks=checks,improved_endpoints=improved,strong_checks=r['checks'])
 result=dict(candidate='interval25',strong_support=bool(strong),experimental_release_allowed=bool(safe),experimental=not bool(strong),broad_strong_checks=broad['checks'],broad_safety_checks=broad_safety,seeds=seeds,known_best_V66_overall=.15088218754984492,predicted_leaderboard_score=None,automatic_upload=False,decision='Proceed with sole candidate as a bounded experimental test; disclose failed strong checks.' if safe else 'Candidate fails prospective experimental safety checks; do not package this candidate.',script_sha256=sha(__file__));write(m.A/'release_decision.json',result);print(result['decision'],flush=True);print('strong support',strong,'experimental allowed',safe,flush=True)
if __name__=='__main__':main()
