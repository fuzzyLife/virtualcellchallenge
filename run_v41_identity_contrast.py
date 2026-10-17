"""Fixed V41 intervention-identity contrast, no new fitted public outcomes."""
import argparse,json,shutil
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v41_identity_contrast as m
import run_v38_standard_counts as counts
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder,aggregate
from audit_transfer_validation import pds,profile
counts.m=m
def prepare():
 m.check();d=m.data();dest=m.A/'predictions';dest.mkdir(exist_ok=True);trace=[]
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);file=dest/(ctx['family']+'.npz');assert not file.exists()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:query=z['query'];measured=z['measured'];base=z['v24'].copy()
  ts=[c['target'] for c in ctx['cells']];ids=[int(np.flatnonzero(d['targets']==t)[0]) for t in ts];base=base[ids];delta,r=m.calibration(key(ctx['family']));effect=m.predict(delta,base,r['use_calibrator'])
  np.savez_compressed(file,v24=base,query=query,measured=measured,identity_contrast25=effect);m.write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'case_targets':ts,**r});trace.append({'family':ctx['family'],**r})
 m.write(m.A/'confidence_gates.json',{'records':trace,'query_outcomes_used_for_abstention':False});print(json.dumps(trace,indent=2),flush=True)
def means():
 counts.means();d=m.data();rows=[];mask=~np.isin(d['genes'],d['targets'])
 for ctx in contexts():
  folder,_=checked_folder(ctx['family'])
  with np.load(m.A/'predictions'/(ctx['family']+'.npz'),allow_pickle=False) as z:effects={a:z[a] for a in ('v24',)+m.ARMS}
  with np.load(folder/'effects.npz',allow_pickle=False) as z:bulk=z['bulk']
  control=profile(sparse.load_npz(folder/'reference.npz'));truth=np.array([profile(sparse.load_npz(folder/f'{i}_actual.npz')) for i in range(len(ctx['cells']))])
  for arm,eff in effects.items():
   pp=[]
   for i,c in enumerate(ctx['cells']):
    p=bulk*np.exp2(eff[i]);p[d['genes']==c['target']]*=.1;p/=p.sum();pp.append(np.log1p(p*50000))
   score=pds(np.array(pp)[:,mask],truth[:,mask],control[mask]);rows.append({'family':ctx['family'],'study':key(ctx['family']),'arm':arm,'pds':score})
 result=m.read(m.A/'mean_results.json');result.update({'pooled_identity_records':rows,'pooled_identity_summary':{a:aggregate(rows,a,'pds') for a in ('v24',)+m.ARMS}});m.write(m.A/'mean_results.json',result);print(json.dumps(result['pooled_identity_summary'],indent=2),flush=True)
def reuse():
 m.check();parent=m.ROOT/'experiments/v38_response_hurdle/analysis';source=m.read(parent/'null_results.json');dest=m.A/'null';dest.mkdir(exist_ok=True);records=[];arm=m.ARMS[0]
 for rr in source['records']:
  if rr['arm'].endswith('_transport'):continue
  aa='v24' if rr['arm']=='v24' else arm;stem=f"{rr['family']}__-1__0__";file=parent/'null'/(stem+rr['arm']+'.npz');assert m.b.sha256(file)==rr['arrays_sha256'];shutil.copyfile(file,dest/(stem+aa+'.npz'));new={**rr,'arm':aa};m.write(dest/(stem+aa+'.json'),new);records.append(new)
 m.write(m.A/'null_results.json',{'records':records,'panels':[],'arms':['v24',arm],'scope':'Exact same zero-effect generator/split/seed, positive weights map zero to zero.','parent_results_sha256':m.b.sha256(parent/'null_results.json')})
 dest=m.A/'development';dest.mkdir(exist_ok=True);reused=[]
 for file in (parent/'development').glob('*__v24.json'):
  rr=m.read(file);assert m.b.sha256(file.with_suffix('.npz'))==rr['arrays_sha256'];shutil.copyfile(file,dest/file.name);shutil.copyfile(file.with_suffix('.npz'),dest/file.with_suffix('.npz').name)
  gate=m.read(m.A/'predictions'/(rr['family']+'.json'))
  if not gate['use_calibrator']:
   stem=file.name.removesuffix('v24.json');new={**rr,'arm':arm};shutil.copyfile(file.with_suffix('.npz'),dest/(stem+arm+'.npz'));m.write(dest/(stem+arm+'.json'),new);reused.append({'family':rr['family'],'case_index':rr['case_index'],'reason':'Source-validation abstention gives bitwise baseline effects and identical emitter/seed.'})
 m.write(m.A/'identical_reuse.json',{'abstaining_treatment_cases':reused,'response_changed_populations_copied':0,'script_sha256':m.b.sha256(m.ROOT/'run_v41_identity_contrast.py')})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','means','null','development','select','confirmation']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='null':reuse()
  elif a.stage in ['development','confirmation']:counts.counts(a.stage)
  elif a.stage=='select':counts.select()
  else:globals()[a.stage]()

