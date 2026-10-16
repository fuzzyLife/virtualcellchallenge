"""Frozen small ensemble experiment with selection before additional count cases."""
import argparse,gc,json,hashlib,platform
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v35_methods as m
from v35_contrast_methods import transform
from v37_ensemble import combine,emit,ARMS
from run_v35_v36 import contexts,key,clean
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics

ROOT=m.ROOT;OUT=ROOT/'experiments/v37_component_ensemble';A=OUT/'analysis';CACHE=A/'cache'
METRICS=('pds','fidelity','reach','jaccard','nmae','common_expression_ratio')

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,v):m.write(p,clean(v))
def digest(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def check():
 p=read(A/'plan.json')
 for name,h in p['dependencies'].items():assert m.sha(ROOT/name)==h,name
 for name,h in p['inputs'].items():assert m.sha(ROOT/name)==h,name
 return p

def freeze():
 A.mkdir(parents=True,exist_ok=True);CACHE.mkdir(exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v37_ensemble.py','run_v37_ensemble.py','test_v37_ensemble.py','v35_methods.py','v36_methods.py','v35_contrast_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py','run_v35_v36.py']
 inputs=['experiments/v34_residual_flow/analysis/plan.json','experiments/v35_response_retrieval/analysis/plan.json','experiments/v24_cd4_response/analysis/leaderboard_result.json','experiments/v36_population_views/analysis/leaderboard_result.json']
 for outer in m.OUTERS:
  for folder in (m.A,m.P36):
   for ext in ('npz','json'):inputs.append(str((folder/'models'/f'{outer}.{ext}').relative_to(ROOT)))
 write(A/'plan.json',{'dependencies':{s:m.sha(ROOT/s) for s in deps},'inputs':{s:m.sha(ROOT/s) for s in inputs},'arms':list(ARMS),'development':'First4 cases/context, original seed. Same public diagnostics previously spent.','confirmation':'Remaining24 raw treatment cases, two seeds (original and +370000). These count-level comparisons are not used to choose the blend; public source and mean labels were previously used. Not new lineage or blind data.','selection':'Require mean MSE/MAE<=1.002 V24; expression ratio<=V24; PDS>=V24-.005; FID>=V24-.002; REACH>=V24; JAC>=V24-.0005; NMAE<=1.002*V24; at least3 of6 metrics improve; every study resolved. Among eligible maximize average oriented relative gain; no retuning. Confirm selected arm with same guards, and per-study PDS/FID/REACH losses <=.02. Otherwise retain V24 and produce no renamed duplicate.','null':'One zero-response control-null comparison/context for each arm; null DE calls cannot exceed V24. Treatment-independent tests only.','context_count':len(contexts()),'case_count':sum(len(c['cells']) for c in contexts()),'D3_opened':False,'python':platform.python_version(),'outcomes_known':'Published V24/V36 and all earlier experiments known. New fixed hypothesis test, no independent private test claimed.'})
 print('V37 plan frozen',flush=True)

def prepare():
 check();m.setup();m.check();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets']
 for ctx in contexts():
  name=ctx['family'];folder=CACHE/name;folder.mkdir(exist_ok=True)
  if (folder/'receipt.json').exists():continue
  outer=key(name);x,meta,measured,source=m.raw(name,genes);controls=x[ctx['left']];depth,cp,bulk,_=moments(controls);query=cp*1e6;query[~measured]=np.nan;base=m.base_for(query,[outer]);net,g,rec=m.load(outer,35)
  pred,_,trace=m.predict(net,g,query,[outer],base=base);shuffle,_,_=m.predict(net,g,query,[outer],base=base,shuffle=True)
  retrieval=transform(pred,base,shuffle,genes,targets,g['axis'],trace,query,'retrieval_contrast25');center=transform(pred,base,shuffle,genes,targets,g['axis'],trace,query,'target_centered25')
  netp,gp,recp=m.load(outer,36);population,_,_=m.predict(netp,gp,query,[outer],base=base,control=m.control_stats(controls));effects=combine(base,retrieval,population,center)
  np.savez_compressed(folder/'effects.npz',**effects,retrieval=retrieval,query=query,cp=cp,bulk=bulk,measured=measured)
  sparse.save_npz(folder/'reference.npz',x[ctx['right']]);sparse.save_npz(folder/'null_template.npz',x[ctx['left'][:400]])
  for i,cell in enumerate(ctx['cells']):
   sparse.save_npz(folder/f'{i}_template.npz',x[cell['template']]);sparse.save_npz(folder/f'{i}_actual.npz',x[cell['actual']])
  files={p.name:m.sha(p) for p in folder.glob('*.npz')}
  write(folder/'receipt.json',{'family':name,'excluded_entire_family':outer,'v35_model_sha256':rec['sha256'],'v36_model_sha256':recp['sha256'],'source':source,'files':files,'query_treatments_used':False,'own_targets_zero':bool(all(np.all(effects[a][i,genes==t]==0) for a in effects for i,t in enumerate(targets))),'source_control_rows':ctx['left'],'reference_rows':ctx['right'],'cases':ctx['cells']})
  assert read(folder/'receipt.json')['own_targets_zero'];del x,controls,net,netp;gc.collect();print('V37 cached',name,flush=True)

def checked_folder(name):
 folder=CACHE/name;r=read(folder/'receipt.json')
 for f,h in r['files'].items():assert m.sha(folder/f)==h,(name,f)
 return folder,r

def mean():
 check();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];records=[]
 for ctx in contexts():
  folder,r=checked_folder(ctx['family']);ref=moments(sparse.load_npz(folder/'reference.npz'))[1]*1e6
  with np.load(folder/'effects.npz',allow_pickle=False) as z:
   query=z['query'];measured=z['measured'];effects={a:z[a] for a in ('v24',)+ARMS}
  for i,cell in enumerate(ctx['cells']):
   ti=int(np.flatnonzero(targets==cell['target'])[0]);own=int(np.flatnonzero(genes==cell['target'])[0]);truthmean=moments(sparse.load_npz(folder/f'{i}_actual.npz'))[1]*1e6;valid=measured&(ref>5);valid[own]=False;truth=np.clip(np.log2((truthmean+.1)/(ref+.1)),-3,3);scores={}
   for arm in ('v24',)+ARMS:
    pm=np.nan_to_num(query)*np.exp2(effects[arm][ti]);pm[own]*=.1;pm*=1e6/pm.sum();e=np.clip(np.log2((pm+.1)/(ref+.1)),-3,3)[valid]-truth[valid];scores[arm]={'mse':float(e@e),'mae':float(abs(e).sum()),'cpm':float(np.sum((pm[valid]-truthmean[valid])**2))}
   for arm in ARMS:records.append({'family':ctx['family'],'study':key(ctx['family']),'target':cell['target'],'arm':arm,**{k+'_ratio':scores[arm][k]/scores['v24'][k] for k in scores[arm]}})
 summary={a:{k:aggregate(records,a,k) for k in ('mse_ratio','mae_ratio','cpm_ratio')} for a in ARMS};write(A/'mean_results.json',{'records':records,'summary':summary});print(json.dumps(summary,indent=2),flush=True)

def aggregate(rows,arm,metric):
 study=[]
 for outer in m.OUTERS:
  # Average cases within contexts, then contexts within studies.
  contexts0=sorted({r['family'] for r in rows if r['study']==outer and r['arm']==arm});values=[]
  for family in contexts0:
   vals=[r[metric] for r in rows if r['arm']==arm and r['family']==family and isinstance(r.get(metric),(int,float)) and np.isfinite(r[metric])]
   if vals:values.append(float(np.mean(vals)))
  if not values:return None
  study.append(float(np.mean(values)))
 return float(np.mean(study))

def count(stage):
 check();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];from evaluate_atlas_transfer import load_atlas_model
 atlas,_=load_atlas_model();dest=A/stage;dest.mkdir(exist_ok=True)
 choice=read(A/'selection.json')['selected'] if stage=='confirmation' else None
 if stage=='confirmation' and choice=='v24':print('No eligible ensemble; confirmation not opened',flush=True);return
 arms=('v24',choice) if stage=='confirmation' else ('v24',)+ARMS
 for ctx in contexts():
  name=ctx['family'];folder,r=checked_folder(name);reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:effects={a:z[a].copy() for a in arms};retrieval=z['retrieval'].copy();base=z['v24'].copy();cp=z['cp'].copy();bulk=z['bulk'].copy()
  ids=range(4) if stage=='development' else range(4,len(ctx['cells'])) if stage=='confirmation' else [-1]
  seeds=[0,370000] if stage=='confirmation' else [0]
  for offset in seeds:
   outputs={a:[] for a in arms}
   for i in ids:
    cell=ctx['cells'][i] if i>=0 else ctx['cells'][0];t=cell['target'];ti=int(np.flatnonzero(targets==t)[0]);own=int(np.flatnonzero(genes==t)[0]);template=sparse.load_npz(folder/(f'{i}_template.npz' if i>=0 else 'null_template.npz'));actual=sparse.load_npz(folder/f'{i}_actual.npz') if i>=0 else template;libs=template.sum(1).A1.astype(np.int64);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
    for arm in arms:
     stem=f'{name}__{i}__{offset}__{arm}';file=dest/(stem+'.json');arrays=dest/(stem+'.npz')
     if file.exists():
      rr=read(file);assert m.sha(arrays)==rr['arrays_sha256']
      with np.load(arrays,allow_pickle=False) as zz:profile=zz['pred'].copy();jack=zz['jack'].copy()
     else:
      eff=effects[arm][ti] if i>=0 else np.zeros(len(genes),np.float32);b=base[ti] if i>=0 else eff;ret=retrieval[ti] if i>=0 else eff
      count,diag=emit(template,eff,b,ret,own if i>=0 else None,cell['seed']+ti+offset,(libs,cp,bulk),atlas,arm)
      np.testing.assert_array_equal(count.sum(1).A1,libs);assert count.dtype.kind in 'iu' and count.data.min(initial=0)>=0
      dd=de(count,reference,genes);profile,jack=profile_jackknife_chunked(count);metrics=direction_metrics(dd,truth,genes,own)
      np.savez_compressed(arrays,pred=profile,jack=jack,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
      rr={'family':name,'study':key(name),'target':t,'case_index':i,'seed_offset':offset,'arm':arm,'cells':len(libs),'row_libraries_exact':True,'null_DE_calls':int(np.sum(dd['eligible']&(dd['q']<.05))),**diag,**metrics,'arrays_sha256':m.sha(arrays)};write(file,rr);del count;gc.collect()
     outputs[arm].append((profile,jack,tp,tj,rr));print('V37',stage,name,i,offset,arm,flush=True)
   if stage!='null':
    for arm,oo in outputs.items():
     p,j,tp,tj=(np.array([o[k] for o in oo]) for k in range(4));rr=[o[4] for o in oo];panel=panel_metrics(p,tp,cc,j,tj,cj,genes,[r['target'] for r in rr],targets)
     for metric in ('fidelity','reach','jaccard','nmae'):
      vals=[r[metric] for r in rr if isinstance(r.get(metric),(int,float)) and np.isfinite(r[metric])];panel[metric]=float(np.mean(vals)) if vals else None
     write(dest/f'{name}__{offset}__{arm}__panel.json',{'family':name,'study':key(name),'seed_offset':offset,'arm':arm,**panel})
 result={'panels':[read(p) for p in sorted(dest.glob('*__panel.json'))],'records':[read(p) for p in sorted(dest.glob('*.json')) if not p.name.endswith('__panel.json')],'arms':list(arms),'scope':'Public unscaled diagnostics; no official anchors.'};write(A/(stage+'_results.json'),result)
 print('V37 completed',stage,flush=True)

def guards(summary,base,means,arm,null):
 s=summary[arm];b=base;resolved=all(s[k] is not None and b[k] is not None for k in METRICS)
 if not resolved:return {'resolved':False},False,None
 tests={'resolved':True,'mean_mse':means[arm]['mse_ratio']<=1.002,'mean_mae':means[arm]['mae_ratio']<=1.002,'expression':s['common_expression_ratio']<=b['common_expression_ratio'],'pds':s['pds']>=b['pds']-.005,'fid':s['fidelity']>=b['fidelity']-.002,'reach':s['reach']>=b['reach'],'jac':s['jaccard']>=b['jaccard']-.0005,'nmae':s['nmae']<=1.002*b['nmae'],'null_no_worse':null[arm]}
 better=[s[k]>b[k] if k not in ('nmae','common_expression_ratio') else s[k]<b[k] for k in METRICS];tests['three_metrics_improve']=sum(better)>=3
 gain=float(np.mean([(s[k]-b[k])/max(abs(b[k]),.01)*(1 if k not in ('nmae','common_expression_ratio') else -1) for k in METRICS]))
 return tests,all(tests.values()),gain

def decide(stage):
 check();results=read(A/(stage+'_results.json'));means=read(A/'mean_results.json')['summary'];nulls=read(A/'null_results.json')['records'];basecalls={r['family']:r['null_DE_calls'] for r in nulls if r['arm']=='v24'};null={a:all(r['null_DE_calls']<=basecalls[r['family']] for r in nulls if r['arm']==a) for a in ARMS}
 summaries={a:{k:aggregate(results['panels'],a,k) for k in METRICS} for a in results['arms']};checks={};eligible=[];scores={}
 for arm in results['arms']:
  if arm=='v24':continue
  tests,passed,gain=guards(summaries,summaries['v24'],means,arm,null);checks[arm]=tests;scores[arm]=gain
  if stage=='confirmation':
   tests['worst_study_guard']=all(np.mean([r[k] for r in results['panels'] if r['study']==o and r['arm']==arm])>=np.mean([r[k] for r in results['panels'] if r['study']==o and r['arm']=='v24'])-.02 for o in m.OUTERS for k in ('pds','fidelity','reach'));passed=all(tests.values())
  if passed:eligible.append(arm)
 selected=max(eligible,key=lambda a:scores[a]) if eligible else 'v24';path=A/('selection.json' if stage=='development' else 'confirmation_decision.json');assert not path.exists();write(path,{'selected':selected,'eligible':eligible,'checks':checks,'summaries':summaries,'gain_scores':scores,'results_sha256':m.sha(A/(stage+'_results.json')),'mean_results_sha256':m.sha(A/'mean_results.json'),'null_results_sha256':m.sha(A/'null_results.json'),'confirmed_better_leaderboard':False,'D3_opened':False});print(json.dumps(read(path),indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','mean','null','development','confirmation','select','confirm']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage in ('select','confirm'):decide('development' if a.stage=='select' else 'confirmation')
  elif a.stage in ('null','development','confirmation'):count(a.stage)
  else:globals()[a.stage]()
