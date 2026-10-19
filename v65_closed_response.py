"""V65 available-axis compositional source response transport, anchored to V24 geometry."""
from pathlib import Path
from functools import lru_cache
import argparse,gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v62_conserved_response as old
from v63_identity_guard import guard
from v57_pie import sha,read,write
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v65_closed_response';A=OUT/'analysis'
ARMS=('reclosed50','closed50','closed75')

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v65_closed_response.py','run_v65_experiment.py','test_v65_response.py','v62_conserved_response.py','v63_identity_guard.py','v15_ensemble_methods.py','v7_response_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 write(A/'plan.json',dict(version=65,arms=ARMS,dependencies={n:sha(ROOT/n) for n in deps},source_cache_sha256=sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),history_sha256=sha(A/'published_history.json'),hypothesis='V63 improves expression but loses fidelity/reach. Test source available-gene-axis normalization and within-source perturbation closure before log1p assay conversion. Compare reclosure alone, closed composition at half correction, and closed composition at three-quarter correction; all retain V24 margin guard and the unchanged residual10 emitter.',estimator='Close each finite source CPM control profile over available genes to 1e6. For closed arms form source perturbed profile with source LFC and own-target factor 0.1, then close it on that same available axis. Transform log1p(0.05*perturbed_CPM)-log1p(0.05*control_CPM), subtract leave-one-target-out source-common templates (>=20 other supported targets), use inherited uncertainty weighting/shrinkage and query inversion. Missing source effects assume unchanged solely for normalization, never become response evidence. Derivative uncertainty is an approximate diagonal delta method ignoring compositional covariance. This does not recover missing genes or claim full-transcriptome normalization.',selection='Before counts, select lowest equal-study target-weighted desired-profile MSE ratio versus V63 among arms with ratio<1, each study<=1.03, LFC MSE/MAE <=1.02 and equal-study mean PDS no loss>.01. Ties by ARMS order. Then lock; compare generated cells with V63, V24 and selected source-target cyclic derangement. No post-count candidate promotion. Complete all fixed unique-target cases at two seeds as diagnostics, not unseen biology. A package may be explicitly experimental on public count guard failure but must state it. No guaranteed private improvement.',count_guards=dict(profile_MSE_ratio=1.0,pds_loss=.005,fidelity_loss=.002,reach_loss=.005,jaccard_loss=.0005,nmae_ratio=1.01),validation='Canonical first case per target/context from existing V34/V37 pools: 30 unique target/context pairs, five contexts, three studies. Generation seeds 0 and620000 are technical repeats. Public biology reused; D3 unopened. Complete held-study source exclusion (HCT116 and HEK293T jointly). No query treatment in fitting.',D3_opened=False,automatic_upload=False))
 print('V65 frozen',sha(A/'plan.json'),flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==p['source_cache_sha256']
 assert sha(A/'published_history.json')==p['history_sha256'];return p

def source_arrays(bank,compositional):
 q=np.maximum(np.nan_to_num(bank['profiles'],nan=0.),0.);q=q/np.maximum(q.sum(1,keepdims=True),1e-30)*1e6
 factor=np.exp2(np.where(bank['ok'],bank['y'],0.))
 if compositional:factor[:,np.arange(len(bank['own'])),bank['own']]=.1
 changed=q[:,None,:]*factor
 if compositional:changed*=1e6/np.maximum(changed.sum(2,keepdims=True),1e-30)
 assay=np.log1p(.05*changed)-np.log1p(.05*q[:,None,:])
 shared=(bank['ok'].sum(2)>0).sum(0)>=2;mask=bank['ok']&shared[None,:,None]
 n=mask.sum(1)[:,None,:]-mask;total=np.where(mask,assay,0).sum(1)[:,None,:]
 template=np.divide(total-np.where(mask,assay,0),n,out=np.zeros_like(assay),where=n>=20)
 derivative=np.log(2)*(.05*changed)/(1+.05*changed)
 variance=np.maximum(bank['var']*derivative**2,.0025)
 return assay-template,variance,q

@lru_cache(8)
def pooled(excluded,shuffle=False):
 b=old.bank(excluded);idx=np.roll(np.arange(300),1) if shuffle else np.arange(300);dest={}
 for kind in ('reclosed','closed'):
  y,v,q=source_arrays(b,kind=='closed')
  dest[kind]=old.pool(y[:,idx],v[:,idx],b['ok'][:,idx],old.WEIGHTS,b['own'])
  del y,v;gc.collect()
 return dest

def predict(query,baseline,bulk,excluded=(),shuffle=False):
 b=old.bank(tuple(sorted(excluded)));pooled_effect=pooled(tuple(sorted(excluded)),shuffle);own=b['own'];panel=np.isin(old.data()['genes'],old.data()['targets']);idx=np.roll(np.arange(300),1) if shuffle else np.arange(300)
 supported=b['ok'][:,idx].any(0)&np.isfinite(query)[None,:];out={};traces={}
 for arm in ARMS:
  dest=old.invert_assay(pooled_effect['reclosed' if arm=='reclosed50' else 'closed'],query)
  dest=np.asarray([old.gated(e.copy(),np.nan_to_num(query,nan=0.),int(o)) for e,o in zip(dest,own)])
  strength=.75 if arm=='closed75' else .5
  proposed=np.where(supported,np.clip(baseline+strength*np.clip(dest-baseline,-1,1),-3,3),baseline).astype(np.float32);proposed[np.arange(300),own]=0
  out[arm],traces[arm]=guard(baseline,proposed,bulk,own,panel)
 return out,traces

def prepare():
 check();from run_v35_v36 import contexts,key
 from run_v37_ensemble import checked_folder
 for ctx in contexts():
  name=ctx['family'];path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  folder,_=checked_folder(name)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();q=z['query'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
  pred,traces=predict(q,base,bulk,(key(name),));sh,shtrace=predict(q,base,bulk,(key(name),),True)
  vp=ROOT/'experiments/v63_identity_guard/analysis/predictions'/(name+'.npz')
  assert sha(vp)==read(vp.with_suffix('.json'))['sha256']
  with np.load(vp,allow_pickle=False) as z:v63=z['all_effects'].copy()
  np.savez_compressed(path,query=q,bulk=bulk,measured=measured,v24=base,v63=v63,**pred,**{k+'_shuffled':v for k,v in sh.items()})
  write(path.with_suffix('.json'),dict(sha256=sha(path),traces=traces,shuffled_traces=shtrace,source_receipt_sha256=sha(folder/'receipt.json'),v63_prediction_sha256=sha(vp),excluded_study=key(name)))
  old.bank.cache_clear();pooled.cache_clear();gc.collect();print('V65 prepared',name,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);a=p.parse_args()
 with threadpool_limits(limits=2):globals()[a.stage]()
