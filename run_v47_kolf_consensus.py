"""Frozen KOLF corroboration validation on identical V24 public splits."""
import argparse,shutil
import numpy as np
from threadpoolctl import threadpool_limits
import v47_kolf_consensus as m
import run_v47_count_grid as counts
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
counts.m=m
def prepare():
 m.check();d=m.data();dest=m.A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);file=dest/(ctx['family']+'.npz');assert not file.exists()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();query=z['query'];measured=z['measured']
  ts=[c['target'] for c in ctx['cells']];ids=[int(np.flatnonzero(d['targets']==t)[0]) for t in ts];pred=m.predict(base);np.savez_compressed(file,v24=base[ids],query=query,measured=measured,**{a:v[ids] for a,v in pred.items()})
  m.write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'case_targets':ts,'held_source_study':key(ctx['family']),'independent_KOLF_study_allowed':True,'changed_outputs':{a:int((pred[a]!=base).sum()) for a in pred}})
def reuse():
 m.check();src=m.ROOT/'experiments/v38_response_hurdle/analysis';dest=m.A/'null';dest.mkdir(exist_ok=True);records=[]
 for rr in m.read(src/'null_results.json')['records']:
  if rr['arm'].endswith('_transport'):continue
  arms=['v24'] if rr['arm']=='v24' else m.ARMS;file=src/'null'/f"{rr['family']}__-1__0__{rr['arm']}.npz";assert m.b.sha256(file)==rr['arrays_sha256']
  for a in arms:
   stem=f"{rr['family']}__-1__0__{a}";shutil.copyfile(file,dest/(stem+'.npz'));new={**rr,'arm':a};m.write(dest/(stem+'.json'),new);records.append(new)
 m.write(m.A/'null_results.json',{'records':records,'panels':[],'arms':['v24',*m.ARMS],'zero_response_exact_V24':True,'parent_sha256':m.b.sha256(src/'null_results.json')})
 dest=m.A/'development';dest.mkdir(exist_ok=True)
 for file in (src/'development').glob('*__v24.json'):
  rr=m.read(file);assert m.b.sha256(file.with_suffix('.npz'))==rr['arrays_sha256'];shutil.copyfile(file,dest/file.name);shutil.copyfile(file.with_suffix('.npz'),dest/file.with_suffix('.npz').name)
def select():
 counts.select();s=m.read(m.A/'selection.json');a=m.ARMS[0];tests=s['checks'][a];tests['mean_gate']=a in m.read(m.A/'mean_results.json')['eligible'];tests['beats_matched_shuffle']=s['relative_gain'][a]>s['relative_gain'][m.ARMS[1]]
 base=s['summary']['v24'];p=s['summary'][a];tests['two_metrics_improve']=sum(p[k]<base[k] if k in ('nmae','common_expression_ratio') else p[k]>base[k] for k in counts.METRICS)>=2
 s['eligible']=[a] if all(tests.values()) else [];s['selected']=a if s['eligible'] else 'v24';m.write(m.A/'selection.json',s);print('V47 selection',s['selected'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','means','null','select']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='freeze':m.freeze()
  elif args.stage=='null':reuse()
  elif args.stage=='means':counts.means()
  else:globals()[args.stage]()

