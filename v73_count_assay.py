"""V73: null-centered count likelihood and bulk-assay response transfer."""
from functools import lru_cache
import argparse,gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v72_batch_offset as base
import v70_count_teacher as old
from v57_pie import read,write,sha
ROOT=base.ROOT;OUT=ROOT/'experiments/v73_count_assay';A=OUT/'analysis'
SOURCES=base.SOURCES;ARMS=('raw25','assay25','assay_center25','assay_center50','assay_credible25')
data=base.data;contexts=base.contexts;cases=base.cases;study=base.study;checked_folder=base.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v73_count_assay.py','v72_batch_offset.py','v71_evidence_projection.py','v70_count_teacher.py','v16_methods.py','run_v70_counts.py']
 write(A/'plan.json',dict(version=73,dependencies={n:sha(ROOT/n) for n in deps},sources=SOURCES,arms=ARMS,hypothesis='V72 restores source coverage but raw proportional transfer fails. Use a unit-response-centered log-normal prior on the count likelihood, then transfer source responses on the actual bulk log1p(50000*p) assay scale that motivated successful V63/V66. Retain raw-scale control and shuffled evidence.',estimator='Same batch-matched expected count E, observed K, empirical overdispersion Fano and propagated NTC uncertainty. Optimize log response eta under quasi-Poisson(K;E*exp(eta))/Fano with Normal(0,(2*ln2)^2) prior; 20 bounded Newton steps. At K=E the posterior mode is exactly zero even with weak evidence. Curvature supplies working uncertainty; no exact Bayesian coverage claim. Source gene support and batch rules unchanged from V72.',assay='Save matched baseline proportion E/sum(E) per source/target. Transform source LFC to delta log1p(50000*p), propagate variance by derivative, optionally subtract leave-one-target-out common response with >=20 other targets. Cap precision per study, shrink to zero, include between-study heterogeneity for uncertainty projection. Invert onto query BULK control proportions, apply fixed bounded quarter/half correction to V66, own zero; unchanged residual10 emitter.',selection='Same first-three distinct cases/context and >=0.5% profile gain, no study >2% worse, genuine better than matched cyclic source shuffle. Development repeatedly reused. Lock selected arm before remaining cases. No failed-guard fallback.',confirmation='Same unchanged paired count evaluator and broader >=30-cell diagnostic. Production requires generated-profile/worst-study/shuffle and broad gates. Other failures and unavailable metrics explicit; no hidden score forecast.',V72_result_sha256=sha(base.A/'development_mean_results.json'),D3_opened=False,automatic_upload=False))
 print('V73 null-centered count/assay protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 return p

def fit_log_response(k,e,fano):
 eta=np.zeros_like(k,dtype=float);precision=1/(2*np.log(2))**2
 for _ in range(20):
  mu=e*np.exp(eta);gradient=(k-mu)/fano-precision*eta;curvature=mu/fano+precision
  eta=np.clip(eta+gradient/curvature,-3*np.log(2),3*np.log(2))
 variance=1/(e*np.exp(eta)/fano+precision)
 return eta/np.log(2),variance/np.log(2)**2

def response(x,ids,groups,control):
 k=np.asarray(x[ids].sum(0)).ravel();e=np.zeros(x.shape[1]);cvtotal=np.zeros_like(e);cross=np.zeros_like(e);square=np.zeros_like(e)
 for batch in np.unique(groups[ids]):
  selected=ids[groups[ids]==batch];xx=x[selected];depth=np.asarray(xx.sum(1)).ravel();cp,cv=control[batch];exposure=depth.sum()
  e+=exposure*cp;cvtotal+=exposure**2*cv;cross+=cp*np.asarray(xx.T@depth).ravel();square+=cp*cp*(depth@depth)
 gain=np.divide(k,e,out=np.zeros_like(k),where=e>0)
 residual=np.maximum(np.asarray(x[ids].power(2).sum(0)).ravel()-2*gain*cross+gain*gain*square,0)
 n=len(ids);fano=np.clip(residual*n/(n-1)/np.maximum(k,1),1,1000)
 effect,var=fit_log_response(k,e,fano)
 var+=np.divide(cvtotal,e*e,out=np.full_like(e,1e6),where=e>0)/np.log(2)**2
 return effect,var,k,e

def fit_source(name):
 check();d=data();genes=d['genes'];targets=d['targets'];path=A/'sources'/(name+'.npz');path.parent.mkdir(exist_ok=True)
 if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];return
 from v25_methods import raw
 x,meta,measured,receipt=raw(name,genes);x=sparse.csr_matrix(x,dtype=float);label=meta.gene_target.to_numpy(str);batch=meta['sample'].astype(str).to_numpy() if 'sample' in meta else np.full(len(meta),'unknown');ids=np.flatnonzero(label=='Non-Targeting');control={}
 for b in np.unique(batch[ids]):
  selected=ids[batch[ids]==b]
  if len(selected)>=20:control[b]=base.control_rate(x[selected])
 shape=(300,len(genes));y=np.zeros(shape,np.float32);v=np.full(shape,1e6,np.float32);ok=np.zeros(shape,bool);bulk=np.zeros(shape,np.float32);rows=[]
 for ti,target in enumerate(targets):
  ids=np.flatnonzero(label==target);ids=ids[np.isin(batch[ids],list(control))]
  if len(ids)<40:continue
  effect,var,k,e=response(x,ids,batch,control);mask=measured&(k>=20)&(e>=20);own=list(genes).index(target);mask[own]=False;effect[own]=0
  y[ti]=effect;v[ti]=var;ok[ti]=mask;bulk[ti]=e/e.sum()
  rows.append(dict(target=str(target),cells=len(ids),batches=len(np.unique(batch[ids])),supported_genes=int(mask.sum())))
 np.savez_compressed(path,effect=y,variance=v,valid=ok,baseline_bulk=bulk,genes=genes,targets=targets)
 write(path.with_suffix('.json'),dict(sha256=sha(path),source=name,study=study(name),targets=rows,source_record=receipt,plan_sha256=sha(A/'plan.json')))
 print('V73 source',name,len(rows),flush=True);del x;gc.collect()

@lru_cache(9)
def teacher(excluded,representation):
 shape=(300,18533);groups={}
 for name in SOURCES:
  if study(name) in excluded:continue
  path=A/'sources'/(name+'.npz');r=read(path.with_suffix('.json'));assert sha(path)==r['sha256']
  with np.load(path,allow_pickle=False) as z:y=z['effect'].astype(float);v=z['variance'].astype(float);ok=z['valid'].copy();bulk=z['baseline_bulk'].astype(float)
  if representation!='raw':
   term=50000*bulk*np.exp2(y);derivative=np.log(2)*term/(1+term)
   y=np.log1p(term)-np.log1p(50000*bulk);v=np.maximum(v*derivative**2,.0025)
  if representation=='center':
   n=ok.sum(0)[None,:]-ok;total=np.where(ok,y,0).sum(0)
   template=np.divide(total[None,:]-np.where(ok,y,0),n,out=np.zeros_like(y),where=n>=20)
   y-=template;ok=ok&(n>=20)
  w=ok/(v+.25**2)
  if study(name) not in groups:groups[study(name)]=[np.zeros(shape),np.zeros(shape)]
  groups[study(name)][0]+=w*y;groups[study(name)][1]+=w
 num=np.zeros(shape);den=np.zeros(shape);parts=[]
 for n,w in groups.values():
  mean=np.divide(n,w,out=np.zeros(shape),where=w>0);weight=np.minimum(w,16);num+=weight*mean;den+=weight;parts.append((mean,weight))
 raw=np.divide(num,den,out=np.zeros(shape),where=den>0);mean=num/(den+4);variance=1/(den+4)
 variance+=np.divide(sum(w*(y-raw)**2 for y,w in parts),den,out=np.zeros(shape),where=den>0)
 confidence=(den>0)*np.minimum(1,np.sqrt(den)/4)
 return mean.astype(np.float32),variance.astype(np.float32),confidence.astype(np.float32)

def invert_assay(delta,bulk):
 changed=np.expm1(np.clip(np.log1p(50000*bulk)[None,:]+delta,0,25))/50000
 ratio=np.divide(changed,bulk[None,:],out=np.ones_like(changed),where=bulk[None,:]>0)
 return np.clip(np.log2(np.maximum(ratio,1/8)),-3,3)

def predict(parent,bulk,excluded=(),shuffle=False):
 d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);out={}
 for arm in ARMS:
  representation='raw' if arm=='raw25' else ('center' if arm.startswith('assay_center') else 'assay')
  mean,var,weight=teacher(tuple(sorted(excluded)),representation)
  if shuffle:mean=np.roll(mean,1,axis=0);var=np.roll(var,1,axis=0);weight=np.roll(weight,1,axis=0)
  destination=mean if representation=='raw' else invert_assay(mean,bulk)
  if arm=='assay_credible25':
   width=1.96*np.sqrt(var);lo=invert_assay(mean-width,bulk);hi=invert_assay(mean+width,bulk);destination=np.clip(parent,lo,hi)
  alpha=.5 if arm=='assay_center50' else .25
  e=np.clip(parent+alpha*weight*np.clip(destination-parent,-.75,.75),-3,3).astype(np.float32);e[np.arange(300),own]=0;e[:,bulk==0]=parent[:,bulk==0]
  assert np.max(abs(e-parent))<=alpha*.75+1e-6;out[arm]=e
 return out

def prepare_all():
 check()
 for ctx in contexts():
  name=ctx['family'];path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  parent=ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
  with np.load(parent,allow_pickle=False) as z:b=z['assay75'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
  pred=predict(b,bulk,(study(name),));sh=predict(b,bulk,(study(name),),True)
  np.savez(path,v66=b,bulk=bulk,query=query,measured=measured,**pred,**{a+'_shuffled':e for a,e in sh.items()})
  write(path.with_suffix('.json'),dict(sha256=sha(path),excluded_study=study(name),parent_sha256=sha(parent),plan_sha256=sha(A/'plan.json')))
  print('V73 prepared',name,flush=True);gc.collect()

def means(stage):
 oa,ar,ck=old.A,old.ARMS,old.check
 try:old.A=A;old.ARMS=ARMS;old.check=check;old.mean_screen(stage)
 finally:old.A=oa;old.ARMS=ar;old.check=ck

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','prepare','development','confirmation']);p.add_argument('--source');args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='fit':fit_source(args.source)
  elif args.stage=='prepare':prepare_all()
  elif args.stage in ('development','confirmation'):means(args.stage)
  else:freeze()
