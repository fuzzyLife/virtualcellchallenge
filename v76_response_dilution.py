"""V76: composition-aware response dilution including direct target knockdown."""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as old
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v76_response_dilution';A=OUT/'analysis'
ARMS=('mixture075','mixture050','dose030','dose055')
data=old.data;contexts=old.contexts;cases=old.cases;study=old.study;checked_folder=old.checked_folder

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v76_response_dilution.py','v16_methods.py','v70_count_teacher.py','run_v75_verified_counts.py']
 inputs={f'experiments/v66_assay_strength/analysis/predictions/{c["family"]}.npz':sha(ROOT/'experiments/v66_assay_strength/analysis/predictions'/(c['family']+'.npz')) for c in contexts()}
 write(A/'plan.json',dict(version=76,dependencies={n:sha(ROOT/n) for n in deps},inputs=inputs,arms=ARMS,hypothesis='V75 improves means but degrades broad PDS: shrinking downstream responses while retaining fixed direct own knockdown changes response composition. Test convex mixtures of the COMPLETE V66 bulk profile and control (response fractions .75/.50), and direct-own-factor .30/.55 controls. No new cellular sampler. Partial-response mixture interpretation is a working biological hypothesis, not identified responder fractions.',representation='Build V66 closed profile with old own factor .1 exactly once. A mixture p=(1-gain)*control+gain*parent is closed and nonnegative before bounds, and jointly dilutes direct knockdown and downstream response. Encode FULL log2(p/control), including the own coordinate. Nonown corrections bounded +/-0.375 versus the closed V66 log-ratio; own allowed to change. All encoded values clipped +/-3, unexpressed query coordinates abstain. Candidate count emission calls the original V16 emitter with own=None so no second suppression is applied. V66 baseline uses its original own argument.',negative_controls='Mixture: cyclic active-target proportion residuals with recipient nonpanel norm; preserve genuine recipient own probability before reclosure. Dose: apply the own-factor change to a cyclic wrong target gene. Decoder constraints can alter exact norm matching; measure final profiles. No outcome-based shuffle choice.',selection='Same first-three distinct development cases/context, >=0.5% profile gain, worst study <=2% worse, beats its shuffled control. Lock before scoring new-model remaining cases. Every public cohort has been used in earlier experiments, including V75 confirmation: no blind claim.',confirmation='PDS requires >=2 distinct targets. Fixed PDS panels use ALL 4/4/6/8/8 unique cases/context (30 total); profile/DE primary summaries use only remaining 1/1/3/5/5 (15 total). This repairs V75 singleton-PDS infeasibility, without dropping the PDS gate. Two seeds 0/+620000; fresh V66 replays and raw-p audits; fixed broad >=30-cell profile/PDS/shuffle gates. Production requires primary/repeat profile/worst-study/PDS/shuffle and broad checks. All other failures/missing remain explicit.',V75_development_sha256=sha(ROOT/'experiments/v75_assay_shrink/analysis/development_mean_results.json'),V75_remaining_mean_sha256=sha(ROOT/'experiments/v75_assay_shrink/analysis/confirmation_mean_results.json'),D3_opened=False,automatic_upload=False))
 print('V76 complete-response dilution frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/n)==h,n
 return p

def parent_prob(effect,bulk,own):
 factor=np.exp2(np.asarray(effect,float));factor[np.arange(len(own)),own]=.1
 p=bulk[None,:]*factor;p/=p.sum(1,keepdims=True)
 return p

def probability(effect,bulk,own,arm,shuffle=False):
 p0=parent_prob(effect,bulk,own);n=len(own);q=bulk[None,:]
 if arm.startswith('mixture'):
  gain=.75 if arm=='mixture075' else .5;genuine=(1-gain)*q+gain*p0
  if not shuffle:return genuine
  delta=p0-q;panel=np.zeros(len(bulk),bool);panel[own]=True;norm=np.linalg.norm(delta[:,~panel],axis=1);active=np.flatnonzero(norm>1e-15)
  if len(active)<2:raise ValueError('Need >=2 active programs for negative control')
  donor=np.roll(active,1);increment=np.zeros_like(delta);increment[active]=delta[donor]*(norm[active]/norm[donor])[:,None]
  p=np.maximum(p0+(gain-1)*increment,0)
  for i in range(n):
   j=own[i];p[i,j]=0;total=p[i].sum();p[i]*=(1-genuine[i,j])/total;p[i,j]=genuine[i,j]
  return p
 factor=np.exp2(np.asarray(effect,float));factor[np.arange(n),own]=.1;target=np.roll(own,1) if shuffle else own
 factor[np.arange(n),target]=.30 if arm=='dose030' else .55
 p=q*factor;p/=p.sum(1,keepdims=True);return p

def predict(parent,bulk,excluded=(),shuffle=False):
 d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);p0=parent_prob(parent,bulk,own)
 baseline=np.log2(np.maximum(np.divide(p0,bulk[None,:],out=np.ones_like(p0),where=bulk[None,:]>0),1e-30));out={}
 for arm in ARMS:
  pp=probability(parent,bulk,own,arm,shuffle)
  e=np.log2(np.maximum(np.divide(pp,bulk[None,:],out=np.ones_like(pp),where=bulk[None,:]>0),1e-30))
  # Bound changes relative to the normalized reference; direct own changes are explicit.
  original_own=e[np.arange(len(own)),own].copy();e=baseline+np.clip(e-baseline,-.375,.375);e[np.arange(len(own)),own]=original_own
  e=np.clip(e,-3,3).astype(np.float32);e[:,bulk==0]=0
  assert np.isfinite(e).all();out[arm]=e
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
  write(path.with_suffix('.json'),dict(sha256=sha(path),parent_sha256=sha(parent),excluded_study=study(name),plan_sha256=sha(A/'plan.json')))
  print('V76 prepared',name,flush=True);gc.collect()

def means(stage):
 from scipy import sparse
 from audit_transfer_validation import profile
 check();d=data();genes=d['genes'];rows=[]
 if stage=='development':arms=('v66',)+ARMS+tuple(a+'_shuffled' for a in ARMS)
 else:selected=read(A/'selection.json')['selected'];assert selected;arms=('v66',selected,selected+'_shuffled')
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name)
  with np.load(A/'predictions'/(name+'.npz'),allow_pickle=False) as z:pred={a:z[a].copy() for a in arms};bulk=z['bulk'].copy()
  ids=cases(ctx)[:3] if stage=='development' else cases(ctx)[3:]
  for i in ids:
   target=ctx['cells'][i]['target'];ti=list(d['targets']).index(target);own=list(genes).index(target);truth=profile(sparse.load_npz(folder/f'{i}_actual.npz'));mask=np.arange(len(genes))!=own
   for arm in arms:
    factor=np.exp2(pred[arm][ti].astype(float))
    if arm=='v66':factor[own]=.1
    p=bulk*factor;p/=p.sum();pp=np.log1p(50000*p)
    rows.append(dict(source=name,study=study(name),target=target,arm=arm,profile_error=float(np.sum((pp[mask]-truth[mask])**2))))
  print('V76 means',stage,name,len(ids),flush=True)
 studies={a:{s:sum(r['profile_error'] for r in rows if r['study']==s and r['arm']==a)/sum(r['profile_error'] for r in rows if r['study']==s and r['arm']=='v66') for s in ('H1','Orion','CD4')} for a in arms if a!='v66'}
 summary={a:float(np.mean(list(v.values()))) for a,v in studies.items()};write(A/(stage+'_mean_results.json'),dict(records=rows,study_ratios=studies,summary=summary))
 if stage=='development':
  eligible=[a for a in ARMS if summary[a]<=.995 and max(studies[a].values())<=1.02 and summary[a]<summary[a+'_shuffled']]
  selected=min(eligible,key=lambda a:summary[a]) if eligible else None;assert not (A/'selection.json').exists();write(A/'selection.json',dict(selected=selected,eligible=eligible,development_sha256=sha(A/'development_mean_results.json')))
 print(json.dumps(dict(stage=stage,summary=summary,selected=selected)),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','development','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='prepare':prepare_all()
  elif args.stage in ('development','confirmation'):means(args.stage)
  else:freeze()
