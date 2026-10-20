"""V75: calibrate trusted V66 response strength directly in bulk assay space."""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as old
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v75_assay_shrink';A=OUT/'analysis'
ARMS=('profile085','profile070','profile050');GAINS=(.85,.70,.50)
data=old.data;contexts=old.contexts;cases=old.cases;study=old.study;checked_folder=old.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v75_assay_shrink.py','v16_methods.py','v70_count_teacher.py','run_v70_counts.py']
 inputs={f'experiments/v66_assay_strength/analysis/predictions/{c["family"]}.npz':sha(ROOT/'experiments/v66_assay_strength/analysis/predictions'/(c['family']+'.npz')) for c in contexts()}
 write(A/'plan.json',dict(version=75,dependencies={n:sha(ROOT/n) for n in deps},inputs=inputs,arms=ARMS,gains=GAINS,hypothesis='V70-V73 replacement response teachers fail. Retain V66 direction/program evidence and test fixed conservative reductions of its bulk log1p(50000*p) response, with the existing cell emitter. Different from V64 source/V24 LFC norm-ratio scaling. Motivated by spent development outcomes; no novelty or blind-validation claim.',calibration='Compute closed profile from V66 effects before own suppression; subtract bulk control profile. Multiply residual by 0.85/0.70/0.50. Negative control adds gain-1 times a cyclic donor residual among active targets, rescaled to the recipient nonpanel response norm; zero-response targets stay zero. At least two active targets are needed for a label derangement. Invert expm1; anchor fold factors at target coordinate to retain stored own effect zero. Unexpressed own uses mean log-factor anchor. Clip change from V66 to +/-0.375 LFC, final +/-3. Exact gain-one and zero-response identities. Own 0.1 applied once by unchanged emitter.',selection='Same first-three distinct cases/context, >=0.5% equal-study profile gain, no study >2% worse, beats corresponding shuffled response increment. Choose smallest profile ratio then lock before remaining cases. No failed-gate fallback.',confirmation='Same paired count and broad all >=30-cell tests, same gates. Production only if generated-profile/worst-study/shuffle/PDS and broad gates pass. Other failed/missing metrics remain explicit. Reused public biology; no hidden overall forecast.',V74_result_sha256=sha(ROOT/'experiments/v74_assay_gain/analysis/development_mean_results.json'),D3_opened=False,automatic_upload=False))
 print('V75 assay-response gain protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/n)==h,n
 return p

def calibrate(effects,bulk,own,gain,shuffle=False):
 effects=np.asarray(effects,np.float32);bulk=np.asarray(bulk,float)
 if gain==1:return effects.copy()
 q=bulk[None,:]*np.exp2(effects.astype(float));q/=q.sum(1,keepdims=True)
 control=np.log1p(50000*bulk)[None,:];residual=np.log1p(50000*q)-control
 increment=residual
 if shuffle:
  panel=np.zeros(len(bulk),bool);panel[own]=True
  norm=np.linalg.norm(residual[:,~panel],axis=1);active=np.flatnonzero(norm>1e-12)
  increment=np.zeros_like(residual)
  if len(active)>=2:
   donor=np.roll(active,1);increment[active]=residual[donor]*(norm[active]/norm[donor])[:,None]
  else:raise ValueError('Insufficient active targets for matched label derangement')
 proposed=np.expm1(np.clip(control+residual+(gain-1)*increment,0,25))/50000
 logfactor=np.log2(np.maximum(np.divide(proposed,bulk[None,:],out=np.ones_like(proposed),where=bulk[None,:]>0),1e-30))
 anchor=logfactor[np.arange(len(own)),own].copy();missing=bulk[own]<=0
 anchor[missing]=np.mean((logfactor-effects)[:,bulk>0],axis=1)[missing]
 candidate=logfactor-anchor[:,None]
 candidate=effects+np.clip(candidate-effects,-.375,.375)
 candidate=np.clip(candidate,-3,3).astype(np.float32);candidate[np.arange(len(own)),own]=0;candidate[:,bulk==0]=effects[:,bulk==0]
 zero=np.max(abs(effects),axis=1)==0
 candidate[zero]=effects[zero]
 assert np.max(abs(candidate-effects))<=.375+1e-6
 return candidate

def predict(parent,bulk,excluded=(),shuffle=False):
 d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']])
 return {a:calibrate(parent,bulk,own,g,shuffle) for a,g in zip(ARMS,GAINS)}

def prepare_all():
 check()
 for ctx in contexts():
  name=ctx['family'];path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  parent=ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
  with np.load(parent,allow_pickle=False) as z:b=z['assay75'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
  pred=predict(b,bulk,(study(name),));sh=predict(b,bulk,(study(name),),True)
  np.savez(path,v66=b,bulk=bulk,query=query,measured=measured,**pred,**{a+'_shuffled':e for a,e in sh.items()})
  write(path.with_suffix('.json'),dict(sha256=sha(path),parent_sha256=sha(parent),excluded_study=study(name),plan_sha256=sha(A/'plan.json')))
  print('V75 prepared',name,flush=True);gc.collect()

def means(stage):
 oa,ar,ck=old.A,old.ARMS,old.check
 try:old.A=A;old.ARMS=ARMS;old.check=check;old.mean_screen(stage)
 finally:old.A=oa;old.ARMS=ar;old.check=ck

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','development','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='prepare':prepare_all()
  elif args.stage in ('development','confirmation'):means(args.stage)
  else:freeze()
