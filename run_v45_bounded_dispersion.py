"""Public count validation for V24 mean-preserving dispersion arms."""
import argparse,shutil
import numpy as np
from threadpoolctl import threadpool_limits
import v45_bounded_dispersion as m
import run_v45_count_grid as counts
from run_v35_v36 import contexts
from run_v37_ensemble import checked_folder
counts.m=m
def prepare():
 m.check();d=m.data();dest=m.A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);file=dest/(ctx['family']+'.npz');assert not file.exists()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();query=z['query'].copy();measured=z['measured'].copy()
  ts=[c['target'] for c in ctx['cells']];ix=[int(np.flatnonzero(d['targets']==t)[0]) for t in ts];base=base[ix]
  np.savez_compressed(file,v24=base,query=query,measured=measured,**{a:base for a in m.ARMS});m.write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'case_targets':ts,'source_receipt_sha256':m.b.sha256(folder/'receipt.json'),'all_mean_effects_bitwise_V24':True})
 for stage in ['null','development']:
  dest=m.A/stage;dest.mkdir(exist_ok=True);src=m.ROOT/'experiments/v38_response_hurdle/analysis'/stage
  for file in src.glob('*__v24.json'):
   rr=m.read(file);assert m.b.sha256(file.with_suffix('.npz'))==rr['arrays_sha256'];shutil.copyfile(file,dest/file.name);shutil.copyfile(file.with_suffix('.npz'),dest/file.with_suffix('.npz').name)
 print('V45 exact V24 means prepared',flush=True)
def select():
 counts.select();s=m.read(m.A/'selection.json');eligible=[]
 for arm in m.ARMS:
  base=s['summary']['v24'];p=s['summary'][arm]
  n=sum(p[k]<base[k] if k in ['nmae','common_expression_ratio'] else p[k]>base[k] for k in counts.METRICS);s['checks'][arm]['two_metrics_improve']=n>=2
  if all(s['checks'][arm].values()):eligible.append(arm)
 s['eligible']=eligible;s['selected']=max(eligible,key=lambda a:s['relative_gain'][a]) if eligible else 'v24';m.write(m.A/'selection.json',s);print('V45 selection',s['selected'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','means','null','development','select','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='freeze':m.freeze()
  elif args.stage in ['null','development','confirmation']:counts.counts(args.stage)
  elif args.stage=='means':counts.means()
  else:globals()[args.stage]()


