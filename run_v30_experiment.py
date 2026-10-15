"""Freeze and execute the V30 full-target mean screen before generated-count testing."""
from pathlib import Path
import argparse,json,time,sys,importlib.metadata
import numpy as np
from threadpoolctl import threadpool_limits
from concurrent.futures import ProcessPoolExecutor,as_completed
import v30_methods as m
import build_submission as b
from v25_methods import baseline

DEPENDENCIES=['v30_methods.py','run_v30_experiment.py','test_v30_bilinear.py','fit_v7_response.py','v25_methods.py','v24_methods.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py','v8_methods.py','build_submission.py']
INPUTS=['reference/K562_gwps_raw_bulk_01.h5ad','reference/v13_resources/rpe1_raw_bulk_01.h5ad','reference/v30_structured_response/gene_string_embeddings.v0.3.h5','experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v24_cd4_response/analysis/models/teacher_all.npz']

def freeze():
 m.A.mkdir(parents=True,exist_ok=True);p=m.A/'plan.json';assert not p.exists()
 m.write(p,{'version':30,'dependencies':{n:b.sha256(m.ROOT/n) for n in DEPENDENCIES},'inputs':{n:b.sha256(m.ROOT/n) for n in INPUTS},'ranks':[16,32],'regularization':[1,10],'ridge_penalty':'lambda * source fitting targets /64','blend':[.125,.25],'correction_LFC_cap':.5,'outer_families':m.OUTERS,'grouping':m.FAMILY,'variants':m.VARIANTS,'primary_candidates':['string_full','construct_full'],'selection':'For each outer family, minimize family-balanced LFC squared-error ratio over other available H1/Orion/CD4 inner families, with BOTH outer and inner excluded. Choose rank/lambda/alpha only using inner outcomes. Primary count comparators use matching selected full STRING hyperparameters for panel-only and randomized-embedding controls. Construct candidate independently inner-selected.','mean_gate':'Before any generated-count or reserve work, a primary candidate must improve family-balanced masked LFC MSE by>=1% and masked LFC MAE by>=0.5% versus held-family V24; worst-family MSE ratio<=1.02; beat matched random and panel-only variants on both MSE and MAE. Full true-DE effects are not used to pick magnification.','scope':'Exploratory on already-studied public families; no fresh blind claims. Whole-Orion (HCT116+HEK293T) excluded together. Whole-CD4 excluded across all conditions. Replogle K562/RPE1 are one family vote. Source-specific output and embedding preprocessing never use another family. Own targets excluded from exact bilinear loss; missing outputs abstain.','guide_scope':'Published pseudobulk construct labels and fold_expr available; distinct promoter/transcript constructs are heterogeneous interventions, not individual guide replicates. Construct random-effects mean only when>=2 effective constructs. Single-construct labels unchanged. Query efficiency remains unobserved.','conditional_work':'Count tests and subsequent learnable-mean population transport/full transition adaptation require a passing mean finalist. No donor reserve, packaging or upload before count and reserve gates.','runtime':{'python':sys.version,'versions':{p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','scikit-learn','h5py','anndata','threadpoolctl']}}});print('V30 frozen',b.sha256(p),flush=True)

def fit_job(pair):
 with threadpool_limits(limits=2):m.fit_one(*pair)
 return pair

def train():
 plan=m.load_plan()
 for p,h in plan['inputs'].items():assert b.sha256(m.ROOT/p)==h
 todo=[]
 for name in m.SOURCES:
  for variant in m.VARIANTS:
   f=m.A/'models'/f'{name}_{variant}.npz'
   if f.exists():assert b.sha256(f)==json.loads(f.with_suffix('.json').read_text())['sha256']
   else:todo.append((name,variant))
 with ProcessPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(fit_job,t) for t in todo]):print('Completed',f.result(),flush=True)

def rows_for(data,key):
 values=['HCT116','HEK293T'] if key=='Orion' else [key]
 return np.flatnonzero(np.isin(data['row_families'],values))

def summarize_context(data,rows,pred,base):
 num=den=ma=mb=0.;valid_count=0
 for ri in rows:
  ti=int(data['row_targets'][ri]);valid=data['valid'][ri].copy();own=np.flatnonzero(data['genes']==data['targets'][ti])[0];valid[own]=False
  if valid.sum()<100:continue
  truth=data['truth'][ri,valid];a=pred[ti,valid]-truth;c=base[ti,valid]-truth
  num+=float(a@a);den+=float(c@c);ma+=float(abs(a).sum());mb+=float(abs(c).sum());valid_count+=int(valid.sum())
 assert valid_count and den>0 and mb>0
 return {'mse_ratio':num/den,'mae_ratio':ma/mb,'squared_error':num,'baseline_squared_error':den,'absolute_error':ma,'baseline_absolute_error':mb,'valid_coordinates':valid_count,'rows':len(rows)}

def summaries():
 m.load_plan();d=m.read_cache();all_records=[];bases={};graphs={}
 def evaluate(excluded,test,variant,config):
  r,l,a=config;recs=[]
  for ci in np.unique(d['row_contexts'][rows_for(d,test)]):
   rows=rows_for(d,test);rows=rows[d['row_contexts'][rows]==ci];query=d['profiles'][ci];bk=(tuple(excluded),int(ci))
   if bk not in bases:bases[bk]=baseline(d,query,m.excluded_legacy(excluded))
   base=bases[bk];key=(tuple(excluded),int(ci),variant,r,l)
   if key not in graphs:
    est,mask,detail=m.predict_graph(query,excluded,variant,r,l);own=np.array([np.flatnonzero(d['genes']==t)[0] for t in d['targets']]);est=np.array([m.gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(est,own)]);delta=np.clip(est-base,-.5,.5)*mask;delta[np.arange(300),own]=0;graphs[key]=delta.astype(np.float32)
   pred=np.clip(base+a*graphs[key],-3,3);recs.append(summarize_context(d,rows,pred,base))
  return {'mse_ratio':float(np.mean([x['mse_ratio'] for x in recs])),'mae_ratio':float(np.mean([x['mae_ratio'] for x in recs])),'contexts':recs}
 choices={}
 for outer in m.OUTERS:
  choices[outer]={}
  for variant in ['string_full','construct_full']:
   fits=[]
   for config in m.CONFIGS:
    values=[]
    for inner in m.OUTERS:
     if inner==outer:continue
     res=evaluate([outer,inner],inner,variant,config);all_records.append({'stage':'inner','outer':outer,'test':inner,'variant':variant,'config':list(config),**res});values.append(res['mse_ratio'])
    fits.append((float(np.mean(values)),config))
   choices[outer][variant]=list(min(fits,key=lambda v:(v[0],v[1]))[1]);print('V30 inner selection',outer,variant,choices[outer][variant],flush=True)
  for variant in m.VARIANTS:
   config=choices[outer]['construct_full' if variant=='construct_full' else 'string_full'];res=evaluate([outer],outer,variant,config);all_records.append({'stage':'outer','outer':outer,'test':outer,'variant':variant,'config':config,**res});print('V30 outer',outer,variant,res['mse_ratio'],res['mae_ratio'],flush=True)
  graphs.clear();bases.clear()
 outer=[r for r in all_records if r['stage']=='outer'];summary={v:{'mse_ratio':float(np.mean([r['mse_ratio'] for r in outer if r['variant']==v])),'mae_ratio':float(np.mean([r['mae_ratio'] for r in outer if r['variant']==v])),'worst_family_mse_ratio':max(r['mse_ratio'] for r in outer if r['variant']==v)} for v in m.VARIANTS};passes={}
 for v in ['string_full','construct_full']:
  t=summary[v];matched=[r for r in outer if r['variant']==v]
  # Construct-specific negative controls must use ITS choices, not the STRING winner's.
  if v=='construct_full':
   comparators={}
   for control in ['string_panel','random_full']:
    vals=[evaluate([o],o,control,choices[o][v]) for o in m.OUTERS];comparators[control]={k:float(np.mean([x[k] for x in vals])) for k in ['mse_ratio','mae_ratio']}
  else:comparators={c:summary[c] for c in ['string_panel','random_full']}
  checks={'mse_gain':t['mse_ratio']<=.99,'mae_gain':t['mae_ratio']<=.995,'worst_family':t['worst_family_mse_ratio']<=1.02,'beats_random':all(t[k]<comparators['random_full'][k] for k in ['mse_ratio','mae_ratio']),'beats_panel_only':all(t[k]<comparators['string_panel'][k] for k in ['mse_ratio','mae_ratio'])};passes[v]={'passes':all(checks.values()),'checks':checks,'matched_controls':comparators}
 m.write(m.A/'mean_screen.json',{'plan_sha256':b.sha256(m.A/'plan.json'),'choices':choices,'records':all_records,'summary':summary,'release_mean_gates':passes,'new_mean_candidate_passes':any(v['passes'] for v in passes.values()),'scope':'Public reused development; masked LFC errors are mean-screen quantities, not official anchored leaderboard scores. No generated cells or reserve scored.'});print(json.dumps({'summary':summary,'gates':passes},indent=2),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','train','summary']);a=p.parse_args()
 with threadpool_limits(limits=2):
  {'freeze':freeze,'prepare':m.prepare,'train':train,'summary':summaries}[a.stage]()
if __name__=='__main__':main()
