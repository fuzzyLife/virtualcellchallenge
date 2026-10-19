"""Paired profile-error and seed/influence sensitivity; no confidence claims."""
import argparse,importlib
import numpy as np
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v65_experiment import cases

def score(rows,arm):
 study=[]
 for group in ('H1','Orion','CD4'):
  sources=sorted({r['source'] for r in rows if r['study']==group});means={a:[] for a in ('v63',arm)}
  for source in sources:
   for a in means:means[a].append(float(np.mean([r['error'] for r in rows if r['source']==source and r['arm']==a])))
  study.append(np.mean(means[arm])/np.mean(means['v63']))
 return float(np.mean(study))

def main(version):
 m=importlib.import_module('v65_closed_response' if version==65 else 'v66_assay_strength');r=read(m.A/'count_assessment.json');arm=r['selected'];d=m.old.data();genes=d['genes'];rows=[]
 for ctx in contexts():
  for i in cases(ctx):
   target=ctx['cells'][i]['target'];own=list(genes).index(target)
   for offset in (0,620000):
    for a in ('v63',arm):
     p=m.A/'counts'/f"{ctx['family']}__{i}__{offset}__{a}.json";receipt=read(p);assert sha(p.with_suffix('.npz'))==receipt['arrays_sha256']
     with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:axis=np.arange(len(genes))!=own;err=float(np.sum((z['pred'][axis]-z['truth'][axis])**2))
     rows.append(dict(source=ctx['family'],study=key(ctx['family']),target=target,seed_offset=offset,arm=a,error=err))
 overall=score(rows,arm);assert abs(overall-r['paired_equal_study_profile_error_ratio'])<1e-12
 ids=sorted({(v['source'],v['target']) for v in rows});loo=[dict(source=s,target=t,ratio=score([v for v in rows if (v['source'],v['target'])!=(s,t)],arm)) for s,t in ids]
 seeds={str(seed):score([v for v in rows if v['seed_offset']==seed],arm) for seed in (0,620000)}
 write(m.A/'paired_profile_sensitivity.json',dict(candidate=arm,ratio=overall,seed_ratios=seeds,leave_one_context_target_out=loo,loo_range=[min(v['ratio'] for v in loo),max(v['ratio'] for v in loo)],auditor_sha256=sha(__file__),scope='Exactly the frozen equal-context/equal-study plug-in expression ratio. Technical-seed and leave-one-context-target-out influence checks; no confidence interval, independent biological replication or private anchored score.'))
 print('V'+str(version),'paired profile ratio',overall,'seed ratios',seeds,'LOO range',min(v['ratio'] for v in loo),max(v['ratio'] for v in loo),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--version',type=int,choices=[65,66],required=True);main(p.parse_args().version)
