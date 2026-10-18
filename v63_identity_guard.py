"""V63: retain V24 target margins while accepting useful assay-space correction."""
from pathlib import Path
import argparse,gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v62_conserved_response as parent
from v57_pie import sha,read,write
ROOT=parent.ROOT;OUT=ROOT/'experiments/v63_identity_guard';A=OUT/'analysis';ARMS=('guarded_assay',)


def probabilities(effect,bulk,own):
 f=np.exp2(np.asarray(effect,float));f[np.arange(len(own)),own]=.1
 p=f*bulk[None,:];return p/p.sum(1,keepdims=True)


def guard(base,new,bulk,own,panel):
 p0=probabilities(base,bulk,own);pn=probabilities(new,bulk,own);control=np.log1p(50000*bulk);axis=~panel
 d0=(np.log1p(50000*p0)-control)[:,axis];norm=np.linalg.norm(d0,axis=1);u=np.divide(d0,norm[:,None],out=np.zeros_like(d0),where=norm[:,None]>0)
 s0=d0@u.T;original=np.diag(s0)[:,None]-s0;floor=.8*np.maximum(original,0);np.fill_diagonal(floor,-np.inf)
 alpha=np.ones(len(own));bad=np.zeros(len(own),bool)
 for step in range(14):
  p=p0+alpha[:,None]*(pn-p0);factor=np.divide(p,bulk[None,:],out=np.ones_like(p),where=bulk[None,:]>0)
  scale=np.divide(bulk[own]*.1,p[np.arange(len(own)),own],out=np.ones(len(own)),where=p[np.arange(len(own)),own]>0);factor*=scale[:,None]
  e=np.clip(np.log2(np.maximum(factor,1e-30)),-3,3).astype(np.float32);e[np.arange(len(own)),own]=0;e[:,bulk==0]=base[:,bulk==0];e[alpha==0]=base[alpha==0]
  # Check the actual float32 effects and the closed normalized profile that emission receives.
  pp=probabilities(e,bulk,own);delta=(np.log1p(50000*pp)-control)[:,axis];s=delta@u.T;newmargin=np.diag(s)[:,None]-s
  tolerance=1e-9*np.maximum(norm,1e-8)[:,None];bad=np.any(newmargin+tolerance<floor,axis=1)
  if not bad.any():break
  if step==12:alpha[bad]=0
  else:alpha[bad]*=.5
 assert not bad.any(),'Identity-margin guard did not find feasible interpolation'
 return e,dict(full_correction_targets=int(np.sum(alpha==1)),fallback_targets=int(np.sum(alpha==0)),minimum_fraction=float(alpha.min()),fractions=alpha.tolist(),retained_baseline_margin_fraction=.8,all_target_margin_constraints_passed=True)


def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v63_identity_guard.py','run_v63_counts.py','test_v63_guard.py','v62_conserved_response.py','v62_pool_cache.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 write(A/'plan.json',dict(candidate='guarded_assay',dependencies={n:sha(ROOT/n) for n in deps},parent_plan_sha256=sha(parent.A/'plan.json'),parent_candidate='assay_center50',known_results='V62 mean screen improves all three studies; partial count panels show HCT small-panel PDS loss, HEK/CD4 case fidelity losses. This is a disclosed outcome-informed follow-up, not blind hypothesis discovery.',hypothesis='Interpolate in normalized bulk proportions towards the fixed V62 assay estimate, accepting the largest dyadic fraction whose realized float32 profile retains at least 80% of every V24 target-versus-target cosine-template margin. This protects prior identity, not truth; no query treatment labels or reference profiles enter. All 300 panel genes excluded from anchor geometry. Exact own-target factor is used only in compositional closure and later applied once by the emitter.',parameters=dict(margin_retention=.8,maximum_halvings=12),selection='One fixed primary; no strength grid selected on outcomes. Evaluate all 20 development cases against V24 and the equally guarded source-target shuffled diagnostic; confirmation remaining24 cases, two seeds, only after development support. Raw expression must improve, PDS loss<=.005, FID loss<=.002, REACH loss<=.005, JAC loss<=.0005, NMAE no worse>1%, >=3 metrics improve and positive oriented relative average; beats shuffled on expression and PDS. Report any failure. Explicit experimental build is permitted solely for the user-requested diagnostic artifact and must not be labeled validated improvement.',D3_opened=False,automatic_upload=False))
 print('V63 guard frozen',sha(A/'plan.json'),flush=True)


def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert sha(parent.A/'plan.json')==p['parent_plan_sha256'];return p


def prepare():
 check();from run_v35_v36 import contexts
 from run_v37_ensemble import checked_folder
 d=parent.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);dest=A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family'];path=dest/(name+'.npz')
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  folder,_=checked_folder(name)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();q=z['query'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
  # V62 saved the exact all-target predictions only for count cohorts. Recreate other target anchors without outcomes.
  from v62_pool_cache import install
  install();ex=__import__('run_v35_v36').key(name);new=parent.predict(q,base,(ex,))['assay_center50'];sh=parent.predict(q,base,(ex,),True)['assay_center50']
  out,trace=guard(base,new,bulk,own,panel);shout,shtrace=guard(base,sh,bulk,own,panel)
  ix=np.array([list(d['targets']).index(c['target']) for c in ctx['cells']]);np.savez_compressed(path,query=q,measured=measured,v24=base[ix],guarded_assay=out[ix],guarded_assay_shuffled=shout[ix],all_effects=out,all_baseline=base)
  write(path.with_suffix('.json'),dict(sha256=sha(path),trace=trace,shuffled_trace=shtrace,source_cache_receipt_sha256=sha(folder/'receipt.json'),parent_plan_sha256=sha(parent.A/'plan.json')));print('V63 prepared',name,'full corrections',trace['full_correction_targets'],'fallbacks',trace['fallback_targets'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);a=p.parse_args()
 with threadpool_limits(limits=2):globals()[a.stage]()
