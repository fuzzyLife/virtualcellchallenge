"""V71: test whether frozen identity geometry blocks count-evidence corrections."""
from pathlib import Path
from functools import lru_cache
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as t
from v57_pie import read,write,sha
ROOT=t.ROOT;OUT=ROOT/'experiments/v71_evidence_projection';A=OUT/'analysis'
ARMS=('raw25','center25','center50','credible25')
data=t.data;contexts=t.contexts;cases=t.cases;study=t.study;checked_folder=t.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists();t.check()
 deps=['v71_evidence_projection.py','v70_count_teacher.py','v16_methods.py','run_v70_counts.py']
 inputs={f'experiments/v70_count_teacher/analysis/sources/{s}.npz':sha(t.A/'sources'/(s+'.npz')) for s in t.SOURCES}
 write(A/'plan.json',dict(version=71,dependencies={n:sha(ROOT/n) for n in deps},inputs=inputs,hypothesis='V70 guarded raw25/raw50 fail development and leave many target responses unchanged. Test raw25 before that guard, count-derived leave-one-target-out common-program removal, and projection of V66 into a source-evidence credible interval. Correct means only, unchanged emitter. No assertion that original self-identity is biological ground truth.',arms=ARMS,selection='Same >=0.5% equal-study desired-profile error improvement, no study >2% worse, beats matched cyclic source-target shuffle on first three unique fixed cases/context. Development already reused in V70. Lock before remaining cases; no failed-guard fallback.',centering='Within each raw source and output gene, subtract other eligible target responses only when >=20 other targets have measured evidence. No recipient outcomes or query target enters its own template.',credible_projection='Move only when V66 lies outside mean +/-1.96*sqrt(posterior variance + between-study heterogeneity). This is a working model interval, not guaranteed cross-context coverage. Quarter-strength step and max 0.75 LFC discrepancy.',response_constraints='Per-coordinate correction <=0.1875 (raw25/center25/credible25) or <=0.375 (center50), own zero, output clip +/-3, unavailable genes abstain. Original V66 identity-margin guard deliberately absent; validate against real held-out responses instead.',confirmation='Reuse run_v70_counts unchanged by binding module to V71. Same paired templates/seeds, V66/shuffle, frozen count gates. Require generated-profile and worst-study/shuffle gates before production; report all remaining failures. Broader fixed all >=30-cell source diagnostic also must pass.',D3_opened=False,automatic_upload=False,V70_development_sha256=sha(t.A/'development_mean_results.json')))
 print('V71 frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/n)==h,n
 return p

@lru_cache(6)
def teacher(excluded,centered):
 shape=(300,18533);groups={}
 for name in t.SOURCES:
  if study(name) in excluded:continue
  p=t.A/'sources'/(name+'.npz');r=read(p.with_suffix('.json'));assert sha(p)==r['sha256']
  with np.load(p,allow_pickle=False) as z:y=z['effect'].astype(float);v=z['variance'].astype(float);ok=z['valid'].copy()
  if centered:
   count=ok.sum(0);total=np.where(ok,y,0).sum(0);n=count[None,:]-ok
   common=np.divide(total[None,:]-np.where(ok,y,0),n,out=np.zeros_like(y),where=n>=20)
   y=y-common;ok=ok&(n>=20)
  w=ok/(v+.25**2)
  if study(name) not in groups:groups[study(name)]=[np.zeros(shape),np.zeros(shape)]
  groups[study(name)][0]+=w*y;groups[study(name)][1]+=w
 num=np.zeros(shape);den=np.zeros(shape);parts=[]
 for n,w in groups.values():
  mean=np.divide(n,w,out=np.zeros(shape),where=w>0);weight=np.minimum(w,16);num+=weight*mean;den+=weight;parts.append((mean,weight))
 raw=np.divide(num,den,out=np.zeros(shape),where=den>0)
 mean=num/(den+4);variance=1/(den+4)
 heterogeneity=sum(w*(value-raw)**2 for value,w in parts)
 variance+=np.divide(heterogeneity,den,out=np.zeros(shape),where=den>0)
 confidence=(den>0)*np.minimum(1,np.sqrt(den)/4)
 return mean.astype(np.float32),variance.astype(np.float32),confidence.astype(np.float32)

def predict(base,bulk,excluded=(),shuffle=False):
 out={};d=data();own=np.array([list(d['genes']).index(g) for g in d['targets']])
 for arm in ARMS:
  mean,var,weight=teacher(tuple(sorted(excluded)),arm.startswith('center'))
  if shuffle:mean=np.roll(mean,1,axis=0);var=np.roll(var,1,axis=0);weight=np.roll(weight,1,axis=0)
  if arm=='credible25':
   width=1.96*np.sqrt(var);destination=np.clip(base,mean-width,mean+width)
  else:destination=mean
  alpha=.5 if arm=='center50' else .25
  e=np.clip(base+alpha*weight*np.clip(destination-base,-.75,.75),-3,3).astype(np.float32)
  e[np.arange(300),own]=0;e[:,bulk==0]=base[:,bulk==0]
  assert np.max(abs(e-base))<=alpha*.75+1e-6
  out[arm]=e
 return out

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
  print('V71 prepared',name,trace,flush=True);gc.collect()

def means(stage):
 # Reuse the exact frozen V70 mean evaluation and selection, with a separate module binding.
 oldA,oldarms,oldcheck=t.A,t.ARMS,t.check
 try:t.A=A;t.ARMS=ARMS;t.check=check;t.mean_screen(stage)
 finally:t.A=oldA;t.ARMS=oldarms;t.check=oldcheck

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','development','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='prepare':prepare_all()
  elif args.stage in ('development','confirmation'):means(args.stage)
  else:freeze()
