"""V53: V52 with train-constant numerical and missingness channels disabled.

Preserves the V52 experiment; shares immutable labels and fixed training settings.
"""
from functools import lru_cache
import argparse,copy,json
from pathlib import Path
import numpy as np,torch
from threadpoolctl import threadpool_limits
import v52_joint_response as m
import run_v52_joint_response as run
BASE=m.ROOT/'experiments/v52_joint_mean_rank/analysis'
A=m.ROOT/'experiments/v53_joint_guard/analysis'

def normalizer(x):
 med=np.nan_to_num(np.nanmedian(x,axis=0));v=np.where(np.isfinite(x),x,med);raw=np.std(v,axis=0)
 scale=np.where(raw>1e-8,np.maximum(raw,.05),np.inf)
 return med.astype(np.float32),scale.astype(np.float32)

def transform(x,med,scale):
 missing=~np.isfinite(x);active=np.isfinite(scale);v=np.where(missing,med,x)
 return np.column_stack([np.clip((v-med)/scale,-8,8),missing.astype(np.float32)*active]).astype(np.float32)

@lru_cache(5)
def labels(source):
 p=BASE/'labels'/(source+'.npz');rec=m.read(p.with_suffix('.json'))
 assert rec['plan_sha256']==m.b.sha256(BASE/'plan.json') and m.b.sha256(p)==rec['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def install():
 m.A=A;m.normalizer=normalizer;m.transform=transform;m.labels=labels

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 p=copy.deepcopy(m.read(BASE/'plan.json'));p['version']=53
 p['dependencies'].update({n:m.b.sha256(m.ROOT/n) for n in ['run_v53_joint_guard.py','test_v53_joint_guard.py']})
 p['inputs']['experiments/v52_joint_mean_rank/analysis/plan.json']=m.b.sha256(BASE/'plan.json')
 for source in m.SOURCES:
  path=BASE/'labels'/(source+'.json');p['inputs'][str(path.relative_to(m.ROOT))]=m.b.sha256(path)
 p['correction']='Exact train-constant features (std<=1e-8) and their missingness indicators are fixed to zero at fitting and inference; unsupported donor slots cannot activate untrained weights. Same model, labels, folds, loss, seeds and selection; all 12 outer networks retrained.'
 p['known_predecessor_evidence']=['V52 source-validation histories observed','V52 H1 held-study metrics observed','V52 other held-study population artifacts may already exist; no biological threshold or hyperparameter changed']
 p['predecessor_prediction_receipts_existing_at_freeze']={f.name:m.b.sha256(f) for f in (BASE/'predictions').glob('*.json')}
 p['scope']='Numerical feature-contract correction to V52, selected from implementation review after partial predecessor outcomes. Not a blind new hypothesis or a private leaderboard claim.'
 m.write(A/'plan.json',p);folder=A/'labels';folder.mkdir(exist_ok=True)
 for source in m.SOURCES:
  r=m.read(BASE/'labels'/(source+'.json'));r['parent_label_plan_sha256']=r['plan_sha256'];r['plan_sha256']=m.b.sha256(A/'plan.json');r['linked_npz']=str((BASE/'labels'/(source+'.npz')).relative_to(m.ROOT));m.write(folder/(source+'.json'),r)
 print('V53 frozen',m.b.sha256(A/'plan.json'),flush=True)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','train','evaluate','audit']);parser.add_argument('--outer',choices=[*m.OUTERS,'all']);parser.add_argument('--kind',choices=['labels','boundaries','artifacts'],default='artifacts');args=parser.parse_args()
 install();torch.set_num_threads(2);torch.use_deterministic_algorithms(True)
 with threadpool_limits(limits=2):
  if args.stage=='freeze':freeze()
  elif args.stage=='train':
   for outer in ([args.outer] if args.outer else [*m.OUTERS,'all']):run.train(outer)
  elif args.stage=='evaluate':run.evaluate()
  else:
   import audit_v52_joint_response as audit
   {'labels':audit.label_replay,'boundaries':audit.boundaries,'artifacts':audit.artifacts}[args.kind]()
if __name__=='__main__':main()
