"""Release only prospectively selected and confirmed V25 geometry responses."""
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy.stats import rankdata
from v25_methods import A,ELIGIBLE,load_plan,write
import build_submission as b
from assess_v24_experiment import compare as old_compare

def compare(f,name,guards):
 base=f.copy();base['candidate']=base.candidate.replace({'v24':'v16','force25':'cdf25','force50':'cdf50','bures25':'cdf50_de','shuffled25':'shuffled25','shuffled50':'shuffled50'});key={'force25':'cdf25','force50':'cdf50','bures25':'cdf50_de'}[name]
 if name=='bures25':
  rows=f[f.candidate=='shuffled25'].copy();rows['candidate']='shuffled50_de';base=pd.concat([base,rows],ignore_index=True)
 rec=old_compare(base,key,guards);rec['candidate']=name;row=f[f.candidate==name].set_index(['family','seed'])
 for control in ['amplitude25','source_program25']:
  other=f[f.candidate==control].set_index(['family','seed']).loc[row.index];rec['checks']['beats_'+control]=sum(float((row[k]-other[k]).groupby('family').mean().mean())>0 for k in ['fidelity','reach','jaccard'])>=2
 rec['pass']=all(rec['checks'].values());return rec

def main(stage):
 plan=load_plan();dest=A/('selection.json' if stage=='development' else 'confirmation_result.json');assert not dest.exists();e=json.loads((A/f'{stage}_evaluation.json').read_text());assert e['complete'] and e['panels_sha256']==b.sha256(A/f'{stage}_panels.csv');f=pd.read_csv(A/f'{stage}_panels.csv');names=ELIGIBLE if stage=='development' else [json.loads((A/'selection.json').read_text())['candidate']];rec=[compare(f,n,plan['raw_guards']) for n in names];passing=[v for v in rec if v['pass']];chosen='v24';ranks={}
 if passing:
  keys=['pds','fidelity','reach','jaccard','nmae','relative_expression_worsening'];value=np.array([[v['changes'][k] for k in keys] for v in passing]);value[:,4:]*=-1;rank=np.array([rankdata(-value[:,i]) for i in range(6)]).T.mean(1);ranks={v['candidate']:float(r) for v,r in zip(passing,rank)};chosen=passing[int(np.argmin(rank))]['candidate']
 write(dest,{'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':rec,'metric_improvement_ranks':ranks,'plan_sha256':b.sha256(A/'plan.json'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'scope':'Repeated public biology; new donor reserve fully excludes CD4 source family. Geometry is not causal identification or hidden anchored-score prediction.'});print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
