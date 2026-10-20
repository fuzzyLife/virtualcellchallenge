"""V72: aggregate batch-offset count evidence without small-batch censoring."""
from pathlib import Path
from functools import lru_cache
import argparse,gc,json
import numpy as np
from scipy import sparse
from scipy.special import digamma,polygamma
from threadpoolctl import threadpool_limits
import v70_count_teacher as t
import v71_evidence_projection as prior
from v57_pie import read,write,sha
_SOURCE_TEACHER=prior.teacher.__wrapped__
ROOT=t.ROOT;OUT=ROOT/'experiments/v72_batch_offset';A=OUT/'analysis'
SOURCES=t.SOURCES;ARMS=prior.ARMS
data=t.data;contexts=t.contexts;cases=t.cases;study=t.study;checked_folder=t.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 dependencies=['v72_batch_offset.py','v70_count_teacher.py','v71_evidence_projection.py','v16_methods.py','run_v70_counts.py']
 write(A/'plan.json',dict(version=72,dependencies={n:sha(ROOT/n) for n in dependencies},sources=SOURCES,arms=ARMS,hypothesis='V70/V71 source batch censoring left most Orion targets with no response evidence. Aggregate all eligible treated cells using matched per-batch control-rate/depth offsets before estimating response and uncertainty; do not require ten treated cells in each batch.',estimator='Control Gamma/quasi-Poisson moments per sample with >=20 NTCs. Expected treatment molecules E=sum_b library_exposure_b*control_proportion_b. Aggregate counts K across every matched treatment batch with >=1 cell, requiring >=40 total cells. Estimate overdispersion around depth*control_proportion_b*(K/E). Posterior log response digamma(0.5+K/Fano)-log(E/Fano), variance trigamma plus propagated control-rate variance. Gene support requires >=20 observed and expected molecules. Fano clipped 1..1000; own coordinates excluded. Working approximate count likelihood, not exact independent-cell NB inference.',selection='Same fixed first-three unique public cases/context; >=0.5% equal-study desired-profile error gain, no study >2% worse, better than cyclic label shuffle. Development is reused; lock before remaining cases. No failed-gate fallback.',confirmation='Same paired count evaluator and broader diagnostic as V70/V71. Production requires generated-profile, worst-study, shuffle and broad gates; report all other failures and unavailable values. No change of selected arm after confirmation.',scope='New estimator after failed V70/V71 development. Same public biology; no independent or hidden-test claim. D3 unopened.',V70_result_sha256=sha(t.A/'development_mean_results.json'),V71_result_sha256=sha(prior.A/'development_mean_results.json'),D3_opened=False,automatic_upload=False))
 print('V72 corrected batch-offset estimator frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 return p

def control_rate(x):
 x=sparse.csr_matrix(x,dtype=float);depth=np.asarray(x.sum(1)).ravel();total=np.asarray(x.sum(0)).ravel();exposure=depth.sum();p=total/exposure;n=len(depth)
 residual=np.maximum(np.asarray(x.power(2).sum(0)).ravel()-2*p*np.asarray(x.T@depth).ravel()+p*p*(depth@depth),0)
 fano=np.clip(residual*n/(n-1)/np.maximum(total,1),1,1000)
 shape=.5+total/fano;rate=exposure/fano
 return shape/rate,shape/rate**2

def response(x,ids,groups,control):
 counts=np.asarray(x[ids].sum(0)).ravel();expected=np.zeros(x.shape[1]);control_var=np.zeros_like(expected);cross=np.zeros_like(expected);squared=np.zeros_like(expected)
 for batch in np.unique(groups[ids]):
  selected=ids[groups[ids]==batch];xx=x[selected];depth=np.asarray(xx.sum(1)).ravel();cp,cv=control[batch]
  exposure=depth.sum();expected+=exposure*cp;control_var+=exposure**2*cv
  cross+=cp*np.asarray(xx.T@depth).ravel();squared+=cp*cp*(depth@depth)
 gain=np.divide(counts,expected,out=np.zeros_like(counts),where=expected>0)
 xx2=np.asarray(x[ids].power(2).sum(0)).ravel();residual=np.maximum(xx2-2*gain*cross+gain*gain*squared,0)
 n=len(ids);fano=np.clip(residual*n/(n-1)/np.maximum(counts,1),1,1000)
 shape=.5+counts/fano;effect=(digamma(shape)-np.log(np.maximum(expected,1e-30)/fano))/np.log(2)
 variance=(polygamma(1,shape)+np.divide(control_var,expected*expected,out=np.full_like(expected,1e6),where=expected>0))/np.log(2)**2
 return np.clip(effect,-3,3),variance,counts,expected

def fit_source(name):
 check();d=data();genes=d['genes'];targets=d['targets'];path=A/'sources'/(name+'.npz');path.parent.mkdir(exist_ok=True)
 if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];return
 from v25_methods import raw
 x,meta,measured,receipt=raw(name,genes);x=sparse.csr_matrix(x,dtype=np.float64);label=meta.gene_target.to_numpy(str);batch=meta['sample'].astype(str).to_numpy() if 'sample' in meta else np.full(len(meta),'unknown');ids=np.flatnonzero(label=='Non-Targeting');control={}
 for b in np.unique(batch[ids]):
  selected=ids[batch[ids]==b]
  if len(selected)>=20:control[b]=control_rate(x[selected])
 y=np.zeros((300,len(genes)),np.float32);v=np.full_like(y,1e6);ok=np.zeros_like(y,bool);rows=[]
 for ti,target in enumerate(targets):
  ids=np.flatnonzero(label==target);ids=ids[np.isin(batch[ids],list(control))]
  if len(ids)<40:continue
  effect,var,k,e=response(x,ids,batch,control);mask=measured&(k>=20)&(e>=20);own=list(genes).index(target);mask[own]=False;effect[own]=0
  y[ti]=effect;v[ti]=var;ok[ti]=mask
  rows.append(dict(target=str(target),cells=len(ids),batches=len(np.unique(batch[ids])),supported_genes=int(mask.sum())))
 np.savez_compressed(path,effect=y,variance=v,valid=ok,genes=genes,targets=targets)
 write(path.with_suffix('.json'),dict(sha256=sha(path),source=name,study=study(name),targets=rows,source_record=receipt,plan_sha256=sha(A/'plan.json')))
 print('V72 source',name,'eligible',len(rows),'with evidence',sum(r['supported_genes']>0 for r in rows),flush=True)
 del x;gc.collect()

@lru_cache(6)
def teacher(excluded,centered):
 # Reuse exactly V71's source pooling and centering with only the source-bank binding changed.
 old=prior.t.A
 try:
  prior.t.A=A
  # Calling __wrapped__ avoids contamination from V71's cached arrays.
  return _SOURCE_TEACHER(excluded,centered)
 finally:prior.t.A=old

def predict(base,bulk,excluded=(),shuffle=False):
 old=prior.teacher
 try:prior.teacher=teacher;return prior.predict(base,bulk,excluded,shuffle)
 finally:prior.teacher=old

def prepare_all():
 check()
 for ctx in contexts():
  name=ctx['family'];path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  parent=ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
  with np.load(parent,allow_pickle=False) as z:base=z['assay75'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
  pred=predict(base,bulk,(study(name),));sh=predict(base,bulk,(study(name),),True)
  np.savez(path,v66=base,bulk=bulk,query=query,measured=measured,**pred,**{a+'_shuffled':e for a,e in sh.items()})
  trace={a:dict(changed_targets=int(np.count_nonzero(np.max(abs(e-base),axis=1)>1e-6)),max_change=float(np.max(abs(e-base)))) for a,e in pred.items()}
  write(path.with_suffix('.json'),dict(sha256=sha(path),excluded_study=study(name),parent_sha256=sha(parent),plan_sha256=sha(A/'plan.json'),changes=trace))
  print('V72 prepared',name,trace,flush=True);gc.collect()

def means(stage):
 oldA,oldarms,oldcheck=t.A,t.ARMS,t.check
 try:t.A=A;t.ARMS=ARMS;t.check=check;t.mean_screen(stage)
 finally:t.A=oldA;t.ARMS=oldarms;t.check=oldcheck

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','prepare','development','confirmation']);p.add_argument('--source');args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='fit':fit_source(args.source)
  elif args.stage=='prepare':prepare_all()
  elif args.stage in ('development','confirmation'):means(args.stage)
  else:freeze()
