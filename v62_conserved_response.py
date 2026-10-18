"""V62: distinguish source-common programs from conserved intervention responses."""
from pathlib import Path
from functools import lru_cache
import argparse, hashlib, json, gc, warnings
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from fit_v7_response import read_cache
from v15_ensemble_methods import cube, pool, profiles
from v7_response_methods import FAMILIES
from v8_methods import gated
from v6_methods import moments
from v57_pie import sha, read, write
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v62_conserved_response'; A=OUT/'analysis'
ARMS=('centered50','assay50','assay_center50')
GROUP={'K562':'Replogle','CD4':'CD4','H1':'H1','HCT116':'Orion','HEK293T':'Orion'}
WEIGHTS=np.array([2.,.25,2.,1.,1.])

@lru_cache(1)
def data(): return read_cache()

def freeze():
 A.mkdir(parents=True,exist_ok=True); assert not (A/'plan.json').exists()
 deps=['v62_conserved_response.py','run_v62_counts.py','test_v62_response.py','fit_v7_response.py','v15_ensemble_methods.py','v7_response_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py','v6_methods.py']
 inputs=['experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v7_multisource_adaptive/analysis/source_provenance.json','experiments/v34_residual_flow/analysis/plan.json']
 write(A/'plan.json',dict(version=62,dependencies={n:sha(ROOT/n) for n in deps},inputs={n:sha(ROOT/n) for n in inputs},arms=ARMS,
 hypothesis='Remove donor-common response contamination, or transfer conserved effects in log1p(50000*proportion) assay space rather than proportional LFC space. Independent adaptation motivated by response decomposition; no claim of causal stress identification.',
 estimator='For each excluded outer study, collapse only allowed source responses with existing 90%-dose convention. Source templates use targets supported in at least two allowed source families; missing genes remain absent. Exclude current target from its own template, require at least 20 other valid targets/output. Weighted pooling retains frozen EB and source precision. Assay-space effects use closed source CPM profiles, then log1p(0.05*CPM*2**LFC)-log1p(0.05*CPM); closure occurs on inversion, own factor handled exactly once by unchanged emitter. Half-strength bounded +/-1 LFC correction to fold-specific V24. Three fixed arms, no strength tuning.',
 controls='Selected arm plus cyclic source-target derangement (all eligible source targets shifted by one) in held-study count tests. Centering-only effects use genuine labels; shuffled null is inference ablation, not separately trained model.',
 selection='Choose minimum equal-study raw pooled-profile squared-error ratio among fixed arms with equal-study LFC MSE/MAE <=1.01 and no study pooled-profile error >1.03. If none, report all failures and continue research without claiming improvement. Count development first four cases/context; confirmation remaining cases and second seed +620000 for locked finalist. Require expression improvement, PDS not down >.005, fidelity not down >.002, reach not down >.005, Jaccard not down >.0005, NMAE not worse >1%; at least three of six improve and favorable average. Per-study expression no worse >3%.',
 context_data='Reuse hashed V37 raw populations and V34 fixed case definitions; no newly blind dataset claim. Fold-specific comparator is existing V24 mechanism, not the exact submitted production array. Raw populations 200 or 400 measured cells; never duplicate cells to manufacture power.',
 literature=['https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1','https://github.com/xinyizhanglab/perturbation-decomposition/blob/main/decomposition/anova.py','https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics.md'],D3_opened=False,automatic_upload=False))
 print('V62 frozen',sha(A/'plan.json'),flush=True)

def check():
 p=read(A/'plan.json')
 for k,h in {**p['dependencies'],**p['inputs']}.items(): assert sha(ROOT/k)==h,k
 return p

@lru_cache(8)
def bank(excluded):
 d=data(); ex=[str(f) for f in FAMILIES if GROUP[str(f)] in excluded]
 y,var,ok,own=cube(d,ex); pp=profiles(d); shared=(ok.sum(2)>0).sum(0)>=2
 # Template mask is selected from source support only, never recipient outcomes.
 mask=ok & shared[None,:,None]
 n=mask.sum(1); total=np.where(mask,y,0).sum(1)
 loo_n=n[:,None,:]-mask
 template=np.divide(total[:,None,:]-np.where(mask,y,0),loo_n,out=np.zeros_like(y),where=loo_n>=20)
 centered=y-template
 assay=np.log1p(.05*np.nan_to_num(pp)[:,None,:]*np.exp2(y))-np.log1p(.05*np.nan_to_num(pp)[:,None,:])
 atotal=np.where(mask,assay,0).sum(1)
 atemplate=np.divide(atotal[:,None,:]-np.where(mask,assay,0),loo_n,out=np.zeros_like(assay),where=loo_n>=20)
 return dict(y=y,var=var,ok=ok,own=own,profiles=pp,centered=centered,assay=assay,assay_centered=assay-atemplate,source_templates=template)

def invert_assay(delta,query):
 q=np.nan_to_num(np.asarray(query,float),nan=0.); z=np.log1p(.05*q)[None,:]+delta
 changed=np.expm1(np.clip(z,0,25))/.05
 # Unsupported / unexpressed coordinates abstain; rank routing uses no pseudocount.
 ratio=np.divide(changed,q[None,:],out=np.ones_like(changed),where=q[None,:]>0)
 return np.clip(np.log2(np.maximum(ratio,1/8)), -3,3)

def predict(query,baseline,excluded=(),shuffle=False):
 b=bank(tuple(sorted(excluded))); own=b['own']; ok=b['ok']; var=b['var']; indices=np.roll(np.arange(300),1) if shuffle else np.arange(300)
 out={}; arrays={'centered50':b['centered'],'assay50':b['assay'],'assay_center50':b['assay_centered']}
 for arm,arr in arrays.items():
  if arm=='centered50':
   dest=pool(arr[:,indices],var[:,indices],ok[:,indices],WEIGHTS,own)
  else:
   # Preserve source sampling weights. Pooling in assay space uses converted uncertainties.
   pp=.05*np.nan_to_num(b['profiles'])[:,None,:]
   derivative=np.log(2)*pp*np.exp2(b['y'])/(1+pp*np.exp2(b['y']))
   v=np.maximum(var*derivative**2,.0025)
   dest=invert_assay(pool(arr[:,indices],v[:,indices],ok[:,indices],WEIGHTS,own),query)
  dest=np.asarray([gated(e.copy(),np.nan_to_num(query,nan=0.),int(o)) for e,o in zip(dest,own)])
  supported=ok[:,indices].any(0)&np.isfinite(query)[None,:]
  corrected=np.where(supported,np.clip(baseline+.5*np.clip(dest-baseline,-1,1),-3,3),baseline).astype(np.float32)
  corrected[np.arange(300),own]=0
  out[arm]=corrected
 return out

def prepare():
 check(); from run_v35_v36 import contexts,key
 from run_v37_ensemble import checked_folder
 records=[]; genes=data()['genes'];targets=data()['targets']
 for ctx in contexts():
  name=ctx['family']; folder,receipt=checked_folder(name); output=A/'predictions'/(name+'.npz'); output.parent.mkdir(exist_ok=True)
  if output.exists():assert sha(output)==read(output.with_suffix('.json'))['sha256'];continue
  with np.load(folder/'effects.npz',allow_pickle=False) as z: baseline=z['v24'].copy();query=z['query'].copy();measured=z['measured'].copy()
  pred=predict(query,baseline,(key(name),)); sh=predict(query,baseline,(key(name),),True)
  ix=np.array([list(targets).index(c['target']) for c in ctx['cells']]); np.savez_compressed(output,query=query,measured=measured,v24=baseline[ix],**{k:v[ix] for k,v in pred.items()},**{k+'_shuffled':v[ix] for k,v in sh.items()})
  write(output.with_suffix('.json'),dict(sha256=sha(output),source_cache_receipt_sha256=sha(folder/'receipt.json'),outer=key(name),targets=[c['target'] for c in ctx['cells']],source_families=[str(f) for f in FAMILIES if GROUP[str(f)]!=key(name)]))
  print('Predicted',name,flush=True)

def means():
 check(); from run_v35_v36 import contexts,key
 from run_v37_ensemble import checked_folder
 from audit_transfer_validation import profile
 records=[];genes=data()['genes']
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name);path=A/'predictions'/(name+'.npz'); assert sha(path)==read(path.with_suffix('.json'))['sha256']
  with np.load(path,allow_pickle=False) as z: pred={k:z[k] for k in z.files}
  with np.load(folder/'effects.npz',allow_pickle=False) as z: bulk=z['bulk'].copy()
  q=np.nan_to_num(pred['query']); ref=moments(sparse.load_npz(folder/'reference.npz'))[1]*1e6
  for i,cell in enumerate(ctx['cells']):
   own=list(genes).index(cell['target']); actual=sparse.load_npz(folder/f'{i}_actual.npz');am=moments(actual)[1]*1e6; ap=profile(actual)
   valid=pred['measured']&(ref>5); valid[own]=False;scoreaxis=np.arange(len(genes))!=own
   truth=np.clip(np.log2((am+.1)/(ref+.1)),-3,3);scores={}
   for arm in ('v24',)+ARMS:
    p=q*np.exp2(pred[arm][i]);p[own]*=.1;p*=1e6/p.sum(); delta=np.clip(np.log2((p+.1)/(ref+.1)),-3,3)[valid]-truth[valid]
    # Bulk-query proportion, not cell-mean CPM, matches expression profile statistic.
    bp=bulk*np.exp2(pred[arm][i]);bp[own]*=.1;bp/=bp.sum(); pp=np.log1p(50000*bp)
    scores[arm]=dict(mse=float(delta@delta),mae=float(abs(delta).sum()),profile=float(np.sum((pp[scoreaxis]-ap[scoreaxis])**2)))
   for arm in ARMS:records.append(dict(source=name,study=key(name),target=cell['target'],arm=arm,**{k+'_ratio':scores[arm][k]/scores['v24'][k] for k in scores[arm]}))
  print('Mean screen',name,flush=True)
 from collections import defaultdict
 summary={};studies={}
 for arm in ARMS:
  studies[arm]={s:{k:float(np.mean([r[k] for r in records if r['arm']==arm and r['study']==s])) for k in ('mse_ratio','mae_ratio','profile_ratio')} for s in ('H1','Orion','CD4')}
  summary[arm]={k:float(np.mean([v[k] for v in studies[arm].values()])) for k in ('mse_ratio','mae_ratio','profile_ratio')}
 eligible=[a for a in ARMS if summary[a]['mse_ratio']<=1.01 and summary[a]['mae_ratio']<=1.01 and max(v['profile_ratio'] for v in studies[a].values())<=1.03]
 selected=min(eligible,key=lambda a:summary[a]['profile_ratio']) if eligible else None
 write(A/'mean_results.json',dict(records=records,summary=summary,studies=studies,eligible=eligible,selected=selected,known_reused_public_biology=True))
 print(json.dumps(dict(summary=summary,selected=selected),indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','means']);args=p.parse_args()
 with threadpool_limits(limits=2):globals()[args.stage]()
