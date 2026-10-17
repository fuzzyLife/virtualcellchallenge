"""Independent arithmetic checks for proposed rank-effect labels and external claims."""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json
import numpy as np
from scipy import stats
ROOT=Path(__file__).resolve().parent

def rank_summary(x,c):
    x,c=np.asarray(x,float),np.asarray(c,float)
    pair=x[:,None]-c[None,:]
    u=float((pair>0).sum()+.5*(pair==0).sum())
    _,counts=np.unique(np.r_[x,c],return_counts=True)
    n,m=len(x),len(c);N=n+m
    tie=float(np.sum(counts.astype(float)**3-counts))
    variance=n*m/12*((N+1)-tie/(N*(N-1)))
    p=float(2*stats.norm.sf(abs((u-n*m/2)/np.sqrt(variance)))) if variance>0 else 1.
    return dict(mean=float(np.mean(x)),U=u,rank_effect=u/(n*m)-.5,tie_sum=tie,p=p)

def main():
    fixtures=[([0,0,1,2],[0,0,0,1]),([0,0,0],[0,0,0]),([2,3,4],[0,1,2]),([1,1,1,1],[0,0,2,2])]
    rng=np.random.default_rng(51052)
    fixtures += [(rng.integers(0,9,size=7).tolist(),rng.integers(0,9,size=11).tolist()) for _ in range(30)]
    for x,c in fixtures:
        got=rank_summary(x,c);ref=stats.mannwhitneyu(x,c,method='asymptotic',use_continuity=False)
        assert abs(got['U']-ref.statistic)<1e-12
        assert abs(got['p']-np.nan_to_num(ref.pvalue,nan=1.))<1e-12
    # AUC/U and arithmetic mean do not alone fix tie-corrected p-values.
    control=[0,1,2,3];seen={};counterexample=None
    for x in combinations_with_replacement(range(7),4):
        got=rank_summary(x,control);key=(got['mean'],got['U'])
        if key in seen and abs(got['p']-seen[key][1]['p'])>1e-4:
            counterexample=dict(control=control,population_A=list(seen[key][0]),summary_A=seen[key][1],population_B=list(x),summary_B=got);break
        seen[key]=(x,got)
    assert counterexample is not None
    # Mean direction can oppose rank direction; no decoder may assume equality.
    opposite=rank_summary([0,0,0,20],[1,1,1,1])
    assert opposite['mean']>1 and opposite['rank_effect']<0
    result=dict(passed=True,fixtures_matched_scipy=len(fixtures),same_mean_and_U_different_tie_corrected_p=counterexample,
        mean_up_rank_down_example=opposite,
        conclusion='Predicting mean and rank effect separately is useful; exact p also needs tie structure. Expression score additionally needs pooled profile and jackknife terms. No leaderboard gain established.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out=ROOT/'analysis/v51_plateau/rank_arithmetic_checks.json';out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
