"""Run fingerprint transfer with exact-count comparators and frozen outer exclusions."""
import argparse,shutil,json
import numpy as np
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
from scipy import sparse
import v43_lincs_fingerprints as m
import run_v42_count_grid as counts
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder,aggregate
from audit_transfer_validation import pds,profile
counts.m=m
def job(pair):
 with threadpool_limits(limits=2):m.fit(*pair)
 return pair
def train():
 m.check();m.prepare_sources()
 todo=[]
 for name in m.core.SOURCES:
  for variant in ['string_full','random_full']:
   file=m.A/'models'/(name+'_'+variant+'.npz')
   if file.exists():assert m.b.sha256(file)==m.read(file.with_suffix('.json'))['sha256']
   else:todo.append((name,variant))
 with ProcessPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(job,p) for p in todo]):print('V43 fitted',f.result(),flush=True)
def prepare():
 m.check();d=m.data();dest=m.A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);file=dest/(ctx['family']+'.npz');assert not file.exists()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();query=z['query'].copy();measured=z['measured'].copy()
  ts=[c['target'] for c in ctx['cells']];ix=[int(np.flatnonzero(d['targets']==t)[0]) for t in ts]
  pred,trace=m.predict(query,[key(ctx['family'])],base)
  for arm in pred:pred[arm]=pred[arm][ix]
  np.savez_compressed(file,v24=base[ix],query=query,measured=measured,**pred)
  m.write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'case_targets':ts,'outer':key(ctx['family']),'trace':trace})
  print('V43 predicted',ctx['family'],flush=True)
def means():
 counts.means();d=m.data();records=[];mask=~np.isin(d['genes'],d['targets'])
 for ctx in contexts():
  folder,_=checked_folder(ctx['family'])
  with np.load(m.A/'predictions'/(ctx['family']+'.npz'),allow_pickle=False) as z:effects={a:z[a] for a in ('v24',)+m.ARMS}
  with np.load(folder/'effects.npz',allow_pickle=False) as z:bulk=z['bulk']
  control=profile(sparse.load_npz(folder/'reference.npz'));truth=np.array([profile(sparse.load_npz(folder/f'{i}_actual.npz')) for i in range(len(ctx['cells']))])
  for arm,eff in effects.items():
   pp=[]
   for i,c in enumerate(ctx['cells']):
    p=bulk*np.exp2(eff[i]);p[d['genes']==c['target']]*=.1;p/=p.sum();pp.append(np.log1p(p*50000))
   records.append({'family':ctx['family'],'study':key(ctx['family']),'arm':arm,'pds':pds(np.array(pp)[:,mask],truth[:,mask],control[mask])})
 result=m.read(m.A/'mean_results.json');result.update({'ideal_PDS':{a:aggregate(records,a,'pds') for a in ('v24',)+m.ARMS},'ideal_PDS_records':records});m.write(m.A/'mean_results.json',result);print(json.dumps(result['ideal_PDS'],indent=2),flush=True)
def reuse():
 m.check();src=m.ROOT/'experiments/v38_response_hurdle/analysis';dest=m.A/'null';dest.mkdir(exist_ok=True);records=[]
 for rr in m.read(src/'null_results.json')['records']:
  if rr['arm'].endswith('_transport'):continue
  arms=['v24'] if rr['arm']=='v24' else m.ARMS
  file=src/'null'/f"{rr['family']}__-1__0__{rr['arm']}.npz";assert m.b.sha256(file)==rr['arrays_sha256']
  for arm in arms:
   stem=f"{rr['family']}__-1__0__{arm}";shutil.copyfile(file,dest/(stem+'.npz'));new={**rr,'arm':arm};m.write(dest/(stem+'.json'),new);records.append(new)
 m.write(m.A/'null_results.json',{'records':records,'panels':[],'arms':['v24',*m.ARMS],'scope':'Exact zero-effect standard emitter, unchanged seeds/splits, fingerprint correction not applied to null.','parent_sha256':m.b.sha256(src/'null_results.json')})
 dest=m.A/'development';dest.mkdir(exist_ok=True);reused=[]
 for file in (src/'development').glob('*__v24.json'):
  rr=m.read(file);assert m.b.sha256(file.with_suffix('.npz'))==rr['arrays_sha256'];shutil.copyfile(file,dest/file.name);shutil.copyfile(file.with_suffix('.npz'),dest/file.with_suffix('.npz').name)
  with np.load(m.A/'predictions'/(rr['family']+'.npz'),allow_pickle=False) as z:
   for arm in m.ARMS:
    if not np.array_equal(z[arm][rr['case_index']],z['v24'][rr['case_index']]):continue
    stem=file.name.removesuffix('v24.json');shutil.copyfile(file.with_suffix('.npz'),dest/(stem+arm+'.npz'));m.write(dest/(stem+arm+'.json'),{**rr,'arm':arm});reused.append({'case':stem,'arm':arm})
 m.write(m.A/'identical_reuse.json',{'bitwise_unchanged_cases':reused,'changed_cases_reused':0})
def select():
 counts.select();s=m.read(m.A/'selection.json');mean=m.read(m.A/'mean_results.json');eligible=[]
 for arm in m.ARMS[:2]:
  tests=s['checks'][arm];tests['mean_gate']=arm in mean['eligible'];tests['beats_shuffled_mean']=all(mean['summary'][arm][k]<mean['summary']['shuffled250'][k] for k in ['mse_ratio','mae_ratio'])
  base=s['summary']['v24'];cand=s['summary'][arm];n=sum(cand[k]<base[k] if k in ['nmae','common_expression_ratio'] else cand[k]>base[k] for k in counts.METRICS);tests['two_metrics_improve']=n>=2
  if all(tests.values()):eligible.append(arm)
 s['eligible']=eligible;s['selected']=max(eligible,key=lambda a:s['relative_gain'][a]) if eligible else 'v24';m.write(m.A/'selection.json',s);print('V43 eligible',eligible,'selected',s['selected'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train','prepare','means','null','development','select','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='freeze':m.freeze()
  elif args.stage in ['development','confirmation']:counts.counts(args.stage)
  else:globals()[args.stage]()

