"""V40 identity-focused mean screen followed by unchanged standard count tests."""
import argparse,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v40_identity_calibration as m
import run_v38_standard_counts as count_runner
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder,aggregate
from audit_transfer_validation import pds,profile
count_runner.m=m
def prepare():
 m.check();d=m.data();targets=d['targets'];out=m.A/'predictions';out.mkdir(exist_ok=True)
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);file=out/(ctx['family']+'.npz');assert not file.exists()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:query=z['query'].copy();measured=z['measured'].copy();base=z['v24'].copy()
  ts=[c['target'] for c in ctx['cells']];indices=np.array([np.flatnonzero(targets==t)[0] for t in ts]);base=base[indices]
  theta=m.load(key(ctx['family']));shuffle=m.load(key(ctx['family']),True);effects={a:m.predict(theta,base,float(a.removeprefix('rank'))/100) for a in m.ARMS}
  np.savez_compressed(file,v24=base,query=query,measured=measured,shuffled25=m.predict(shuffle,base,.25),**effects)
  m.write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'case_targets':ts,'outer':key(ctx['family']),'model_sha256':m.b.sha256(m.A/'models'/(key(ctx['family'])+'.npz')),'source_cache_receipt_sha256':m.b.sha256(folder/'receipt.json')})
  print('V40 predicted',ctx['family'],flush=True)
def means():
 count_runner.means();d=m.data();rows=[];genes=d['genes'];nonpanel=~np.isin(genes,d['targets'])
 for ctx in contexts():
  folder,_=checked_folder(ctx['family'])
  with np.load(m.A/'predictions'/(ctx['family']+'.npz'),allow_pickle=False) as z:effects={a:z[a] for a in ('v24',)+m.ARMS+('shuffled25',)}
  with np.load(folder/'effects.npz',allow_pickle=False) as z:bulk=z['bulk']
  control=profile(sparse.load_npz(folder/'reference.npz'));truth=np.array([profile(sparse.load_npz(folder/f'{i}_actual.npz')) for i in range(len(ctx['cells']))])
  for arm,eff in effects.items():
   predictions=[]
   for i,c in enumerate(ctx['cells']):
    prob=bulk*np.exp2(eff[i]);prob[genes==c['target']]*=.1;prob/=prob.sum();predictions.append(np.log1p(prob*50000))
   score=pds(np.array(predictions)[:,nonpanel],truth[:,nonpanel],control[nonpanel]);rows.append({'family':ctx['family'],'study':key(ctx['family']),'arm':arm,'pds':score})
 summary={a:aggregate(rows,a,'pds') for a in ('v24',)+m.ARMS+('shuffled25',)};path=m.A/'mean_results.json';record=m.read(path);eligible=record['eligible']
 selected=max(eligible,key=lambda a:summary[a]) if eligible else max(m.ARMS,key=lambda a:summary[a])
 record.update({'pooled_identity_records':rows,'pooled_identity_summary':summary,'selected':selected,'selection_objective':'Largest pooled-profile PDS among mean MSE/MAE eligible fixed gains. Ideal pooled margins are a proxy, then exact integer count tests.','shuffled_identity_model_is_diagnostic_only':True});m.write(path,record)
 print(json.dumps({'MSE':record['summary'],'pooled_PDS':summary,'selected':selected},indent=2),flush=True)
def reuse_null():
 import shutil
 parent=m.ROOT/'experiments/v38_response_hurdle/analysis';src=m.read(parent/'null_results.json');selected=m.read(m.A/'mean_results.json')['selected'];dest=m.A/'null';dest.mkdir(exist_ok=True);records=[]
 for rr in src['records']:
  if rr['arm'].endswith('_transport'):continue
  arm='v24' if rr['arm']=='v24' else selected;stem=f"{rr['family']}__-1__0__";file=parent/'null'/(stem+rr['arm']+'.npz');assert m.b.sha256(file)==rr['arrays_sha256'];shutil.copyfile(file,dest/(stem+arm+'.npz'));r={**rr,'arm':arm};m.write(dest/(stem+arm+'.json'),r);records.append(r)
 m.write(m.A/'null_results.json',{'records':records,'panels':[],'arms':['v24',selected],'scope':'Exact same zero-effect emitter, split, seed, control fit and template as V38. Positive gene weights map zero to zero. No model-dependent treatments reused.','parent_results_sha256':m.b.sha256(parent/'null_results.json')})
 development=m.A/'development';development.mkdir(exist_ok=True)
 for f in (parent/'development').glob('*__v24.json'):
  r=m.read(f);assert m.b.sha256(f.with_suffix('.npz'))==r['arrays_sha256'];shutil.copyfile(f,development/f.name);shutil.copyfile(f.with_suffix('.npz'),development/f.with_suffix('.npz').name)
 print('V40 zero-effect and V24 comparator reuse checked',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','means','null','development','select','confirmation']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='null':reuse_null()
  elif a.stage in ['development','confirmation']:count_runner.counts(a.stage)
  elif a.stage=='select':count_runner.select()
  else:globals()[a.stage]()

