"""Separate outcome-known corrective experiment; original five arms remain frozen."""
import argparse,gc,json
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v35_methods as m
import run_v37_ensemble as r
from v37_recalibrated import NAME,emit
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked

PARENT=r.A;A=PARENT/'recalibrated';original_folder=r.checked_folder

def check():
 plan=r.read(A/'plan.json')
 for name,h in plan['dependencies'].items():assert m.sha(m.ROOT/name)==h,name
 assert m.sha(PARENT/'plan.json')==plan['parent_plan_sha256']
 return plan

def freeze():
 A.mkdir(exist_ok=True);assert not (A/'plan.json').exists()
 r.write(A/'plan.json',{'parent_plan_sha256':m.sha(PARENT/'plan.json'),'dependencies':{n:m.sha(m.ROOT/n) for n in ['v37_recalibrated.py','run_v37_recalibrated.py','v37_ensemble.py','run_v37_ensemble.py']},'arm':NAME,'known_outcomes':'All44 mean results, all corrected nulls, and the initial H1/HCT116 treatment panels known. Other treatment panels may be available. Correction motivated by row mixture breaking fitted means. No result for this new corrected arm known.','formula':'Same fixed50:50 parent row mixture, then same dual-moment routine and200 iterations using the mean50 response. No extra noise, coefficient or grid.','selection':'Same frozen original mean/null/count thresholds and scoring, no relaxation. Additional count confirmation remains required if selected.','D3_opened':False})
 for ctx in r.contexts():
  src,_=original_folder(ctx['family']);dest=A/'effects'/ctx['family'];dest.mkdir(parents=True)
  with np.load(src/'effects.npz',allow_pickle=False) as z:
   np.savez_compressed(dest/'effects.npz',v24=z['v24'],**{NAME:z['mean50']},retrieval=z['retrieval'],cp=z['cp'],bulk=z['bulk'])
 means=r.read(PARENT/'mean_results.json');r.write(A/'mean_results.json',{'summary':{NAME:means['summary']['mean50']},'source_sha256':m.sha(PARENT/'mean_results.json'),'reason':'Recalibrated mixture targets exactly the half-logFC mean50 moments.'})
 print('Post-mixture calibration experiment frozen',flush=True)

class FolderView:
 def __init__(self,name,original):self.name=name;self.original=original
 def __truediv__(self,name):return A/'effects'/self.name/name if name=='effects.npz' else self.original/name

def folder(name):
 src,rec=original_folder(name);return FolderView(name,src),rec

def install():
 r.A=A;r.check=check;r.checked_folder=folder;r.ARMS=(NAME,);r.emit=emit

def null():
 check();data=m.old.old.read_cache();genes=data['genes'];zero=np.zeros(len(genes),np.float32);from evaluate_atlas_transfer import load_atlas_model
 atlas,_=load_atlas_model();dest=A/'null';dest.mkdir(exist_ok=True);sources=A/'null_corrected';sources.mkdir(exist_ok=True);records=[]
 for ctx in r.contexts():
  name=ctx['family'];n=ctx['null'];files=[dest/f'{name}__-1__0__{arm}.json' for arm in ('v24',NAME)]
  if all(p.exists() for p in files):records.extend(r.read(p) for p in files);continue
  x,meta,measured,source=m.raw(name,genes);template=x[n['template']];actual=x[n['actual']];ref=x[n['reference']];_,cp,bulk,_=moments(x[n['fit_controls']]);libs=template.sum(1).A1.astype(np.int64);truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);cc,cj=profile_jackknife_chunked(ref);r.write(sources/(name+'_source.json'),{'source':source,'split':n})
  for arm in ('v24',NAME):
   count,diag=emit(template,zero,zero,zero,None,n['seed'],(libs,cp,bulk),atlas,arm);dd=de(count,ref,genes);profile,jack=profile_jackknife_chunked(count);metrics=direction_metrics(dd,truth,genes,None);arrays=dest/f'{name}__-1__0__{arm}.npz'
   np.savez_compressed(arrays,pred=profile,jack=jack,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible']);rr={'family':name,'study':r.key(name),'target':'Non-Targeting','case_index':-1,'seed_offset':0,'arm':arm,'cells':len(libs),'row_libraries_exact':bool(np.array_equal(count.sum(1).A1,libs)),'null_DE_calls':int(np.sum(dd['eligible']&(dd['q']<.05))),**diag,**metrics,'arrays_sha256':m.sha(arrays)};assert rr['row_libraries_exact'];r.write(dest/f'{name}__-1__0__{arm}.json',rr);records.append(rr);print('Recalibrated null',name,arm,rr['null_DE_calls'],flush=True);del count;gc.collect()
  del x;gc.collect()
 r.write(A/'null_results.json',{'records':records,'panels':[],'arms':['v24',NAME],'scope':'Same original corrected three-way controls.'})

def development():
 install();r.count('development')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','null','development','select','confirmation','confirm','audit']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  elif a.stage=='null':null()
  elif a.stage=='development':development()
  elif a.stage in ('select','confirm'):install();r.decide('development' if a.stage=='select' else 'confirmation')
  elif a.stage=='confirmation':install();r.count('confirmation')
  else:
   import audit_v37_ensemble as audit
   audit.A=A;audit.check=check;audit.main()
