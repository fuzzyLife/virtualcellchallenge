"""V70: batch-matched raw-count response teacher; unchanged V66 count emitter."""
from pathlib import Path
import argparse,gc,json
import numpy as np
from scipy import sparse
from scipy.special import digamma,polygamma
from threadpoolctl import threadpool_limits
from v57_pie import read,write,sha
from v62_conserved_response import data
from v63_identity_guard import guard
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v70_count_teacher';A=OUT/'analysis'
SOURCES=('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Rest_D2','CD4_Stim48_D2')
ARMS=('raw25','raw50')

def study(name):return key(name)

def freeze():
 A.mkdir(parents=True,exist_ok=True)
 assert not (A/'plan.json').exists()
 dependencies=['v70_count_teacher.py','v16_methods.py','v63_identity_guard.py','v25_methods.py','run_v65_experiment.py','v12_cpu_metrics.py','v8_methods.py']
 write(A/'plan.json',dict(version=70,hypothesis='Improve response means rather than count shape: infer batch-matched bulk proportions directly from raw counts, with depth offsets and moment-matched quasi-Poisson uncertainty. Pool evidence once per study and shrink uncertain responses; correct frozen V66 means at two prespecified strengths. Not an exact NB likelihood or independent biological replication.',sources=SOURCES,arms=ARMS,dependencies={n:sha(ROOT/n) for n in dependencies},protocol=dict(min_treated=40,min_batch_controls=20,min_batch_treated=10,min_total_molecules=20,fano_floor=1,fano_ceiling=1000,gamma_shape_prior=.5,between_study_sd=.25,signal_prior_sd=.5,max_LFC=3,max_correction=.75,identity_margin_fraction=.8),selection='Choose lowest equal-study desired-profile error ratio on first three unique fixed cases/context among raw25/raw50; require <=0.995 and no study >1.02 and genuine evidence better than cyclic target shuffle. If neither passes, report no supported candidate. Lock before remaining unique cases are scored. Reused public studies, not blind validation.',confirmation='Remaining unique cases/context, evaluate desired means then actual integer counts versus V66 and matched shuffled correction using identical control templates and seeds. Profile ratio <=0.995, no study >1.02; PDS decrease <=.005, FID <=.002, REACH <=.005, JAC <=.0005, NMAE no worse 1%. Missing stays missing and fails complete support. Report all failed gates. No unsupported replacement promoted.',production='Only a mean-eligible candidate can be considered; unchanged V66 residual10 generator, original control rows/library sizes and 80% identity margin against V66. Any incomplete count support must remain explicit and experimental.',source_exclusions='Whole H1/Orion/CD4 study excluded from source teacher for recipient evaluation; no recipient treatments/support used by predictor.',D3_opened=False,automatic_upload=False,parent_model_sha256=sha(ROOT/'experiments/v66_assay_strength/analysis/reference_model_arrays.npz')))
 print('V70 frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 return p

def rate(x):
 x=sparse.csr_matrix(x,dtype=float);n=x.shape[0];depth=np.asarray(x.sum(1)).ravel()
 assert n>1 and (depth>0).all()
 total=np.asarray(x.sum(0)).ravel();exposure=depth.sum();p=total/exposure
 # Variance of total counts around a library-size offset; moment-based overdispersion.
 residual=np.asarray(x.power(2).sum(0)).ravel()-2*p*np.asarray(x.T@depth).ravel()+p*p*(depth@depth)
 fano=np.clip(np.maximum(residual,0)*n/(n-1)/np.maximum(total,1),1,1000)
 shape=.5+total/fano
 lograte=digamma(shape)-np.log(exposure/fano)
 variance=polygamma(1,shape)
 return lograte,variance,total,exposure

def fit_source(name):
 check();d=data();genes=d['genes'];targets=d['targets'];path=A/'sources'/(name+'.npz');path.parent.mkdir(exist_ok=True)
 if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];return
 from v25_methods import raw
 x,meta,measured,receipt=raw(name,genes);label=meta.gene_target.to_numpy(str)
 batch=meta['sample'].astype(str).to_numpy() if 'sample' in meta else np.full(len(meta),'unknown')
 controls=np.flatnonzero(label=='Non-Targeting');control={}
 for b in np.unique(batch[controls]):
  ids=controls[batch[controls]==b]
  if len(ids)>=20:control[b]=rate(x[ids])
 y=np.zeros((300,len(genes)),np.float32);var=np.full_like(y,1e6);valid=np.zeros_like(y,bool);rows=[]
 for ti,t in enumerate(targets):
  ids=np.flatnonzero(label==t);ids=ids[np.isin(batch[ids],list(control))]
  if len(ids)<40:continue
  numerator=np.zeros(len(genes));denominator=np.zeros(len(genes));molecules=np.zeros(len(genes));nb=0
  for b in np.unique(batch[ids]):
   selected=ids[batch[ids]==b]
   if len(selected)<10:continue
   tr,tv,tt,_=rate(x[selected]);cr,cv,ct,_=control[b]
   effect=(tr-cr)/np.log(2);noise=(tv+cv)/np.log(2)**2
   mask=measured&(tt>=20)&(ct>=20);weight=np.divide(1,np.maximum(noise,.0025),where=mask,out=np.zeros_like(noise))
   numerator+=weight*np.clip(effect,-3,3);denominator+=weight;molecules+=tt;nb+=1
  mask=denominator>0;yy=np.divide(numerator,denominator,out=np.zeros_like(numerator),where=mask)
  vv=np.divide(1,denominator,out=np.full_like(denominator,1e6),where=mask)
  # A target's own suppression is applied once by the existing emitter.
  own=list(genes).index(t);mask[own]=False;yy[own]=0
  y[ti]=yy;var[ti]=vv;valid[ti]=mask
  rows.append(dict(target=str(t),cells=len(ids),batches=nb,supported_genes=int(mask.sum())))
 np.savez_compressed(path,effect=y,variance=var,valid=valid,genes=genes,targets=targets)
 write(path.with_suffix('.json'),dict(sha256=sha(path),source=name,study=study(name),targets=rows,source_record=receipt,plan_sha256=sha(A/'plan.json')))
 print('V70 source',name,len(rows),flush=True)
 del x;gc.collect()

def pooled(excluded=(),shuffle=False):
 d=data();shape=(300,len(d['genes']));groups={}
 for name in SOURCES:
  if study(name) in excluded:continue
  path=A/'sources'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256']
  with np.load(path,allow_pickle=False) as z:
   np.testing.assert_array_equal(z['genes'],d['genes']);np.testing.assert_array_equal(z['targets'],d['targets'])
   y=z['effect'].astype(float);v=z['variance'].astype(float);ok=z['valid'].copy()
  w=ok/(v+.25**2)
  if study(name) not in groups:groups[study(name)]=[np.zeros(shape),np.zeros(shape)]
  groups[study(name)][0]+=w*y;groups[study(name)][1]+=w
 num=np.zeros(shape);den=np.zeros(shape);votes=np.zeros(shape,int)
 for n,w in groups.values():
  mean=np.divide(n,w,out=np.zeros(shape),where=w>0)
  # A study cannot obtain unlimited precision merely by contributing many conditions.
  weight=np.minimum(w,1/.25**2);num+=weight*mean;den+=weight;votes+=(w>0)
 mean=np.divide(num,den,out=np.zeros(shape),where=den>0)
 variance=np.divide(1,den,out=np.full(shape,1e6),where=den>0)
 mean*=.5**2/(.5**2+variance)
 reliability=(den>0)*np.minimum(1,np.sqrt(den)/4)
 if shuffle:mean=np.roll(mean,1,axis=0);reliability=np.roll(reliability,1,axis=0)
 return mean.astype(np.float32),reliability.astype(np.float32)

def predict(base,bulk,excluded=(),shuffle=False):
 teacher,confidence=pooled(excluded,shuffle);d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);out={}
 for arm,alpha in zip(ARMS,(.25,.5)):
  proposal=np.clip(base+alpha*confidence*np.clip(teacher-base,-.75,.75),-3,3).astype(np.float32);proposal[np.arange(300),own]=0
  out[arm],_=guard(base,proposal,bulk,own,panel)
 return out

def prepare(name):
 check();path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
 if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];return
 parent=ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
 with np.load(parent,allow_pickle=False) as z:base=z['assay75'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
 pred=predict(base,bulk,(study(name),));sh=predict(base,bulk,(study(name),),True)
 np.savez_compressed(path,v66=base,bulk=bulk,query=query,measured=measured,**pred,**{k+'_shuffled':v for k,v in sh.items()})
 write(path.with_suffix('.json'),dict(sha256=sha(path),excluded_study=study(name),parent_sha256=sha(parent),plan_sha256=sha(A/'plan.json')))
 print('V70 prepared',name,flush=True)

def cases(ctx):
 seen=set();out=[]
 for i,c in enumerate(ctx['cells']):
  if c['target'] not in seen:out.append(i);seen.add(c['target'])
 return out

def mean_screen(stage):
 check();from audit_transfer_validation import profile
 d=data();rows=[]
 if stage=='confirmation':selected=read(A/'selection.json')['selected'];assert selected;arms=('v66',selected,selected+'_shuffled')
 else:arms=('v66',)+ARMS+tuple(a+'_shuffled' for a in ARMS)
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name)
  with np.load(A/'predictions'/(name+'.npz'),allow_pickle=False) as z:pred={a:z[a].copy() for a in arms};bulk=z['bulk'].copy()
  ids=cases(ctx);ids=ids[:3] if stage=='development' else ids[3:]
  for i in ids:
   target=ctx['cells'][i]['target'];ti=list(d['targets']).index(target);own=list(d['genes']).index(target);actual=sparse.load_npz(folder/f'{i}_actual.npz');truth=profile(actual);mask=np.arange(len(bulk))!=own
   for arm in arms:
    factor=np.exp2(pred[arm][ti]);factor[own]=.1;p=bulk*factor;p/=p.sum();pp=np.log1p(50000*p)
    rows.append(dict(source=name,study=study(name),target=target,arm=arm,profile_error=float(np.sum((pp[mask]-truth[mask])**2))))
  print('V70 means',stage,name,len(ids),flush=True)
 studies={a:{} for a in arms if a!='v66'}
 for a in studies:
  for s in ('H1','Orion','CD4'):
   rr=[r for r in rows if r['study']==s and r['arm']==a];bb=[r for r in rows if r['study']==s and r['arm']=='v66']
   studies[a][s]=sum(r['profile_error'] for r in rr)/sum(r['profile_error'] for r in bb)
 summary={a:float(np.mean(list(v.values()))) for a,v in studies.items()}
 write(A/(stage+'_mean_results.json'),dict(records=rows,study_ratios=studies,summary=summary))
 if stage=='development':
  eligible=[a for a in ARMS if summary[a]<=.995 and max(studies[a].values())<=1.02 and summary[a]<summary[a+'_shuffled']]
  selected=min(eligible,key=lambda a:summary[a]) if eligible else None
  assert not (A/'selection.json').exists();write(A/'selection.json',dict(selected=selected,eligible=eligible,criterion='Frozen development mean criterion',development_sha256=sha(A/'development_mean_results.json')))
 print(json.dumps(dict(stage=stage,summary=summary,selected=selected)),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','prepare','development','confirmation']);p.add_argument('--source');args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='fit':fit_source(args.source)
  elif args.stage=='prepare':prepare(args.source)
  elif args.stage in ('development','confirmation'):mean_screen(args.stage)
  else:freeze()
