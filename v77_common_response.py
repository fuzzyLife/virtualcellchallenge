"""V77: remove shared predicted trans-response bias, preserve distinct programs."""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as old
import v76_response_dilution as prior
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v77_common_response';A=OUT/'analysis'
ARMS=('common25','common50','common100','projection50')
data=old.data;contexts=old.contexts;cases=old.cases;study=old.study;checked_folder=old.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v77_common_response.py','v16_methods.py','v70_count_teacher.py','v76_response_dilution.py']
 inputs={f'experiments/v66_assay_strength/analysis/predictions/{c["family"]}.npz':sha(ROOT/'experiments/v66_assay_strength/analysis/predictions'/(c['family']+'.npz')) for c in contexts()}
 write(A/'plan.json',dict(version=77,dependencies={n:sha(ROOT/n) for n in deps},inputs=inputs,arms=ARMS,hypothesis='Global dilution trades lower error for worse PDS. Instead shrink only the common predicted trans-response component across the query target panel. Related to V62 source centering, but applied after full V66 prediction in the scoring assay; no source teacher replacement or new sampler.',estimator='LOO average of other predicted log1p(50000*p) responses; exclude every source intervention own coordinate from its output-gene template. Subtract fixed .25/.5/1 times this common component, or half the positively aligned common projection (coefficient clipped0..2). Preserve each target original own probability before nonown mass reclosure. Reanchor with original V66 closure to retain own factor.1, bound changes +/-0.375 LFC and total +/-3, preserve unsupported/zero-trans-response targets. All cells use original emitter with its legacy own argument.',control='Fixed cyclic gene-coordinate permutation of common component within ten bulk-expression strata, separately for panel/nonpanel genes; recipient coefficient reused. Input common-component norms preserved. Final decoder constraints may change achieved norms. No score-based permutation choice.',development='Profile first3 distinct fixed cases/context. PDS first3 plus8 metadata-selected extra targets with >=100 cells, excluding ALL fixed cases from extras, sorted by SHA256(v77-pds:+target). Freeze cohort before computing truth. Joint rank-focused eligibility: PDS gain >=0.01, mean profile ratio <=1.005 (noninferiority), worst study <=1.02, and profile better than gene-shuffled control. Prefer larger PDS gain, tie lower profile error.',validation='Remaining fixed cases for profile/DE; PDS all30 distinct fixed cases to avoid singleton panels; two seeds with fresh V66 replay and raw-p/BH/rank audits. Held-target point check excludes every development-PDS target from all >=30-cell source targets (about600); require PDS gain>=.005, profile ratio<=1.005, worst study<=1.02 and beats gene shuffle. Count production requires PDS gain>=.01, same profile/worst-study/shuffle gates in BOTH seeds. Other DE failures and missing values explicit.',scope='PDS-focused joint experiment with prespecified different objective, not a relaxation/promotion of failed V75/V76 mean candidates. All public biological data were used previously; held targets are new to this coefficient selection, not new biological studies. No hidden scaled-score forecast.',D3_opened=False,automatic_upload=False))
 print('V77 common-response rank-focused protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/n)==h,n
 return p

def permute_common(common,bulk,panel):
 result=common.copy();ix=np.flatnonzero(bulk>0);order=ix[np.argsort(bulk[ix],kind='stable')]
 for group in np.array_split(order,10):
  for ispanel in (False,True):
   subset=group[panel[group]==ispanel]
   if len(subset)>1:result[:,subset]=common[:,np.roll(subset,1)]
 return result

def calibrate(parent,bulk,own,arm,shuffle=False):
 parent=np.asarray(parent,np.float32);bulk=np.asarray(bulk,float);n=len(own);p0=prior.parent_prob(parent,bulk,own)
 response=np.log1p(50000*p0)-np.log1p(50000*bulk)[None,:]
 valid=np.ones(response.shape,bool);valid[np.arange(n),own]=False
 total=np.where(valid,response,0).sum(0);number=valid.sum(0)
 common=(total[None,:]-np.where(valid,response,0))/(number[None,:]-valid)
 panel=np.zeros(len(bulk),bool);panel[own]=True
 if arm=='projection50':
  coefficient=.5*np.clip(np.divide(np.sum(response[:,~panel]*common[:,~panel],axis=1),np.sum(common[:,~panel]**2,axis=1),out=np.zeros(n),where=np.sum(common[:,~panel]**2,axis=1)>0),0,2)
 else:coefficient=np.full(n,{'common25':.25,'common50':.5,'common100':1.}[arm])
 if shuffle:common=permute_common(common,bulk,panel)
 pp=np.expm1(np.clip(np.log1p(50000*p0)-coefficient[:,None]*common,0,25))/50000
 for i,j in enumerate(own):
  pp[i,j]=0;pp[i]*=(1-p0[i,j])/pp[i].sum();pp[i,j]=p0[i,j]
 original_factor=np.exp2(parent.astype(float));original_factor[np.arange(n),own]=.1;closure=np.sum(bulk[None,:]*original_factor,axis=1)
 ratio=np.divide(pp,bulk[None,:],out=np.ones_like(pp),where=bulk[None,:]>0)*closure[:,None]
 proposal=np.log2(np.maximum(ratio,1e-30));e=np.clip(parent+np.clip(proposal-parent,-.375,.375),-3,3).astype(np.float32)
 e[np.arange(n),own]=0;e[:,bulk==0]=parent[:,bulk==0]
 inactive=np.max(abs(parent[:,bulk>0]),axis=1)<=1e-8;e[inactive]=parent[inactive]
 return e

def predict(parent,bulk,excluded=(),shuffle=False):
 d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']])
 return {a:calibrate(parent,bulk,own,a,shuffle) for a in ARMS}

def prepare_all():
 check()
 for ctx in contexts():
  name=ctx['family'];p=A/'predictions'/(name+'.npz');p.parent.mkdir(exist_ok=True)
  if p.exists():assert sha(p)==read(p.with_suffix('.json'))['sha256'];continue
  parent=ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
  with np.load(parent,allow_pickle=False) as z:base=z['assay75'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
  pred=predict(base,bulk,(study(name),));sh=predict(base,bulk,(study(name),),True)
  np.savez(p,v66=base,bulk=bulk,query=query,measured=measured,**pred,**{a+'_shuffled':e for a,e in sh.items()})
  write(p.with_suffix('.json'),dict(sha256=sha(p),parent_sha256=sha(parent),excluded_study=study(name),plan_sha256=sha(A/'plan.json')))
  print('V77 prepared',name,flush=True);gc.collect()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='prepare':prepare_all()
  else:freeze()
