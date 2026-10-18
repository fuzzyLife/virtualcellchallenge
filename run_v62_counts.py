"""V62 frozen paired integer-count evaluation with common controls."""
import argparse,gc,json,time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v62_conserved_response as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder,aggregate
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked
METRICS=('pds','fidelity','reach','jaccard','nmae','common_expression_ratio')

def counts(stage,source=None):
 m.check();choice=read(m.A/'mean_results.json')['selected']; assert choice is not None,'No eligible mean candidate'
 if stage=='confirmation':assert read(m.A/'development_assessment.json')['passed'],'Development did not support confirmation'
 arms=('v24',choice,choice+'_shuffled');atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model()
 d=m.data();genes=d['genes'];targets=d['targets'];dest=m.A/stage;dest.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family']
  if source and source!=name:continue
  folder,receipt=checked_folder(name);reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy();bulk=z['bulk'].copy()
  path=m.A/'predictions'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256']
  with np.load(path,allow_pickle=False) as z:effects={a:z[a].copy() for a in arms}
  ids=range(4) if stage=='development' else range(4,len(ctx['cells']));offsets=(0,) if stage=='development' else (0,620000)
  for offset in offsets:
   outputs={a:[] for a in arms}
   for i in ids:
    cell=ctx['cells'][i];ti=list(targets).index(cell['target']);own=list(genes).index(cell['target']);template=sparse.load_npz(folder/f'{i}_template.npz');actual=sparse.load_npz(folder/f'{i}_actual.npz');libs=template.sum(1).A1.astype(np.int64)
    truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
    for arm in arms:
     file=dest/f'{name}__{i}__{offset}__{arm}.json';arr=file.with_suffix('.npz');begin=time.time()
     if file.exists():
      r=read(file);assert sha(arr)==r['arrays_sha256']
      with np.load(arr,allow_pickle=False) as z:pr=z['pred'].copy();pj=z['jack'].copy()
     else:
      c,diag=emit(template,np.arange(len(libs)),effects[arm][i],own,cell['seed']+ti+offset,'residual10',(libs,cp,bulk),atlas)
      np.testing.assert_array_equal(c.sum(1).A1,libs);assert c.dtype.kind in 'iu' and c.data.min(initial=0)>=0
      dd=de(c,reference,genes);pr,pj=profile_jackknife_chunked(c);metric=direction_metrics(dd,truth,genes,own)
      np.savez_compressed(arr,pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
      r=dict(family=name,study=key(name),target=cell['target'],case_index=i,seed_offset=offset,arm=arm,cells=len(libs),row_libraries_exact=True,**metric,**diag,arrays_sha256=sha(arr),seconds=time.time()-begin);write(file,r);del c;gc.collect()
     outputs[arm].append((pr,pj,tp,tj,r));print('V62',stage,name,i,offset,arm,'FID',round(r['fidelity'],5),'sec',round(time.time()-begin,1),flush=True)
   for arm,vals in outputs.items():
    pr,pj,tp,tj=(np.array([v[k] for v in vals]) for k in range(4));rr=[v[4] for v in vals];p=panel_metrics(pr,tp,cc,pj,tj,cj,genes,[r['target'] for r in rr],targets)
    for k in ('fidelity','reach','jaccard','nmae'):
     vv=[r[k] for r in rr if r.get(k) is not None and np.isfinite(r[k])];p[k]=float(np.mean(vv)) if vv else None
    write(dest/f'{name}__{offset}__{arm}__panel.json',dict(family=name,study=key(name),seed_offset=offset,arm=arm,**p))

def assess(stage):
 choice=read(m.A/'mean_results.json')['selected'];arms=('v24',choice,choice+'_shuffled');dest=m.A/stage
 panels=[read(p) for p in sorted(dest.glob('*__panel.json'))];expected=5*len(arms)*(1 if stage=='development' else 2);assert len(panels)==expected,(len(panels),expected)
 s={a:{k:aggregate(panels,a,k) for k in METRICS} for a in arms};b=s['v24'];v=s[choice];sh=s[choice+'_shuffled']
 finite=all(v[k] is not None and b[k] is not None for k in METRICS)
 tests=dict(resolved=finite,expression=finite and v['common_expression_ratio']<b['common_expression_ratio'],pds=finite and v['pds']>=b['pds']-.005,fid=finite and v['fidelity']>=b['fidelity']-.002,reach=finite and v['reach']>=b['reach']-.005,jac=finite and v['jaccard']>=b['jaccard']-.0005,nmae=finite and v['nmae']<=b['nmae']*1.01)
 if finite:
  improve=[v[k]>b[k] if k not in ('nmae','common_expression_ratio') else v[k]<b[k] for k in METRICS];gain=float(np.mean([(v[k]-b[k])/max(abs(b[k]),.01)*(-1 if k in ('nmae','common_expression_ratio') else 1) for k in METRICS]));tests['three_improve']=sum(improve)>=3;tests['positive_gain']=gain>0
 else: gain=None
 tests['beats_shuffle_expression']=v['common_expression_ratio'] is not None and v['common_expression_ratio']<sh['common_expression_ratio'];tests['beats_shuffle_pds']=v['pds']>sh['pds']
 studies={}
 for study in ('H1','Orion','CD4'):
  vals={a:float(np.mean([p['common_expression_ratio'] for p in panels if p['arm']==a and p['study']==study])) for a in arms};studies[study]=vals
 tests['study_expression_guard']=all(z[choice]<=z['v24']*1.03 for z in studies.values())
 write(m.A/(stage+'_assessment.json'),dict(passed=all(tests.values()),checks=tests,summary=s,studies=studies,relative_gain=gain,selected=choice,panels=panels,D3_opened=False,private_score=None))
 print(json.dumps(dict(passed=all(tests.values()),checks=tests,summary=s),indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation','assess-development','assess-confirmation']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage.startswith('assess-'):assess(a.stage.split('-',1)[1])
  else:counts(a.stage,a.source)
