"""V38 fixed study-excluded tests; raw-count confirmation follows selection."""
import argparse, gc, importlib.metadata, json, platform
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v50_response_prior as m
from run_v35_v36 import contexts, key, clean
from run_v37_ensemble import checked_folder, aggregate
from v6_methods import moments
from v12_cpu_metrics import de, direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
METRICS=('pds','fidelity','reach','jaccard','nmae','common_expression_ratio')
def write(p,v):m.write(p,clean(v))
def check():
 p=m.read(m.A/'plan.json')
 for n,h in p['dependencies'].items():assert m.b.sha256(m.ROOT/n)==h,n
 for n,h in p['inputs'].items():assert m.b.sha256(m.ROOT/n)==h,n
 return p
def freeze():
 m.A.mkdir(parents=True,exist_ok=True);assert not (m.A/'plan.json').exists()
 dep=['v38_response_hurdle.py','run_v38_response_hurdle.py','test_v38_response_hurdle.py','fit_v7_response.py','v15_ensemble_methods.py','v7_response_methods.py','v25_methods.py','v29_methods.py','v29_pairwise_repair.py','v29_swap_kernel.py','v27_methods.py','v26_methods.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py']
 inputs=['experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v7_multisource_adaptive/analysis/source_provenance.json','experiments/v34_residual_flow/analysis/plan.json','experiments/v37_component_ensemble/analysis/ensemble_release.json']
 for c in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+c['family']+'/receipt.json')
 write(m.A/'plan.json',{'dependencies':{s:m.b.sha256(m.ROOT/s) for s in dep},'inputs':{s:m.b.sha256(m.ROOT/s) for s in inputs},'arms':list(m.ARMS),'fixed_parameters':{'pairs_per_recipient_target':256,'hgb_iterations':120,'leaves':15,'class_active_abs_LFC':.25,'class_active_SE_multiplier':2,'weighted_regression_active_weight':8,'hurdle_coefficient':.25,'weighted_coefficient':.25,'selective_coefficient':.5,'selective_probability':.20,'selective_direction_confidence':.8},'selection':'All three mean candidates tested on same spent44 public means. Choose lowest equal-study MSE ratio among candidates with MSE and MAE<=1.01, and test that one plus its empirical-transport emitter on first4 cases/context. Compare all six raw metrics and corrected control-null tests. Confirmation uses remaining24 count cases, original seed and +380000. No claim of blind lineage validation or official anchor replication. No automatic upload.','training':'Whole outer study excluded; additionally exclude whole recipient study from its donor response features. Hash0 targets validation only. Random genes sampled without enrichment, study-balanced weights. K562 additionally used as recipient. Normal/SE labels are diagnostic, not claimed to match competition Wilcoxon.','python':platform.python_version(),'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn','joblib']},'prior_results_known':True,'D3_opened':False})
 print('V38 plan frozen',flush=True)
def prepare():
 check();d=m.data();genes=d['genes'];targets=d['targets'];records=[]
 for ctx in contexts():
  folder,rec=checked_folder(ctx['family']);out=m.A/'predictions';out.mkdir(exist_ok=True);file=out/(ctx['family']+'.npz')
  if file.exists():continue
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();query=z['query'].copy();measured=z['measured'].copy()
  ts=np.array([c['target'] for c in ctx['cells']]);indices=np.array([np.flatnonzero(targets==t)[0] for t in ts])
  pred,trace=m.predict(m.load(key(ctx['family'])),query,base[indices],ts);np.savez_compressed(file,v24=base[indices],query=query,measured=measured,**pred)
  write(file.with_suffix('.json'),{'sha256':m.b.sha256(file),'source_cache_receipt_sha256':m.b.sha256(folder/'receipt.json'),'outer':key(ctx['family']),'case_targets':ts.tolist(),'trace':trace})
  print('V38 predicted',ctx['family'],flush=True)
def means():
 check();d=m.data();genes=d['genes'];records=[]
 for ctx in contexts():
  folder,_=checked_folder(ctx['family']);path=m.A/'predictions'/(ctx['family']+'.npz');assert m.b.sha256(path)==m.read(path.with_suffix('.json'))['sha256']
  with np.load(path,allow_pickle=False) as z:query=z['query'];measured=z['measured'];effects={a:z[a] for a in ('v24',)+m.ARMS}
  ref=moments(sparse.load_npz(folder/'reference.npz'))[1]*1e6
  for i,cell in enumerate(ctx['cells']):
   own=int(np.flatnonzero(genes==cell['target'])[0]);truthmean=moments(sparse.load_npz(folder/f'{i}_actual.npz'))[1]*1e6;valid=measured&(ref>5);valid[own]=False
   truth=np.clip(np.log2((truthmean+.1)/(ref+.1)),-3,3);scores={}
   for arm in ('v24',)+m.ARMS:
    p=np.nan_to_num(query)*np.exp2(effects[arm][i]);p[own]*=.1;p*=1e6/p.sum();delta=np.clip(np.log2((p+.1)/(ref+.1)),-3,3)[valid]-truth[valid]
    scores[arm]={'mse':float(delta@delta),'mae':float(abs(delta).sum()),'cpm':float(np.sum((p[valid]-truthmean[valid])**2))}
   for arm in m.ARMS:records.append({'family':ctx['family'],'study':key(ctx['family']),'target':cell['target'],'arm':arm,**{k+'_ratio':scores[arm][k]/scores['v24'][k] for k in scores[arm]}})
 summary={a:{k:aggregate(records,a,k) for k in ['mse_ratio','mae_ratio','cpm_ratio']} for a in m.ARMS}
 eligible=[a for a in m.ARMS if summary[a]['mse_ratio']<=1.01 and summary[a]['mae_ratio']<=1.01]
 selected=min(eligible,key=lambda a:summary[a]['mse_ratio']) if eligible else min(m.ARMS,key=lambda a:summary[a]['mse_ratio'])
 write(m.A/'mean_results.json',{'records':records,'summary':summary,'eligible':eligible,'selected':selected,'exploratory_if_no_eligible':not bool(eligible)})
 print(json.dumps({'summary':summary,'selected':selected},indent=2),flush=True)
def counts(stage):
 check();d=m.data();genes=d['genes'];targets=d['targets'];from evaluate_atlas_transfer import load_atlas_model
 atlas,_=load_atlas_model();selection=m.read(m.A/'mean_results.json')['selected'];arms=('v24',)+m.ARMS
 if stage=='confirmation':arms=('v24',m.read(m.A/'selection.json')['selected'])
 dest=m.A/stage;dest.mkdir(exist_ok=True);panels=[];records=[]
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name);null=stage=='null'
  if null:
   spec=ctx['null'];fit=set(spec['fit_controls']);right=set(spec['reference']);actualrows=set(spec['actual']);assert not fit&right and not fit&actualrows and not right&actualrows and set(spec['template']).issubset(fit)
   from v25_methods import raw
   x,meta,meas,source=raw(name,genes);reference=x[spec['reference']];template0=x[spec['template']];actual0=x[spec['actual']];_,cp,bulk,_=moments(x[spec['fit_controls']]);del x;gc.collect()
   write(dest/(name+'_source.json'),{'source':source,'split':spec,'actual_reference_fit_disjoint':True})
   ids=[-1];offsets=[0]
  else:
   reference=sparse.load_npz(folder/'reference.npz')
   with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy();bulk=z['bulk'].copy()
   with np.load(m.A/'predictions'/(name+'.npz'),allow_pickle=False) as z:effects={a:z[a].copy() for a in ['v24',*m.ARMS]}
   ids=range(4) if stage=='development' else range(4,len(ctx['cells']));offsets=[0] if stage=='development' else [0,380000]
  cc,cj=profile_jackknife_chunked(reference)
  for offset in offsets:
   oo={a:[] for a in arms}
   for i in ids:
    cell=ctx['cells'][i] if i>=0 else ctx['null'];t=cell['target'];ti=int(np.flatnonzero(targets==t)[0]) if not null else 0;own=int(np.flatnonzero(genes==t)[0]) if not null else None
    template=sparse.load_npz(folder/f'{i}_template.npz') if not null else template0;actual=sparse.load_npz(folder/f'{i}_actual.npz') if not null else actual0;libs=template.sum(1).A1.astype(np.int64);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
    for arm in arms:
     stem=f'{name}__{i}__{offset}__{arm}';file=dest/(stem+'.json');arrays=dest/(stem+'.npz')
     if file.exists():
      rr=m.read(file);assert m.b.sha256(arrays)==rr['arrays_sha256']
      with np.load(arrays,allow_pickle=False) as z:profile=z['pred'].copy();jack=z['jack'].copy()
     else:
      eff=np.zeros(len(genes)) if null else effects[arm][i]
      count,diag=m.emit(template,eff,own,cell['seed']+ti+offset,(libs,cp,bulk),atlas,arm)
      np.testing.assert_array_equal(count.sum(1).A1,libs);assert count.dtype.kind in 'iu' and count.data.min(initial=0)>=0
      dd=de(count,reference,genes);profile,jack=profile_jackknife_chunked(count);metric=direction_metrics(dd,truth,genes,own)
      np.savez_compressed(arrays,pred=profile,jack=jack,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
      rr={'family':name,'study':key(name),'target':t,'case_index':i,'seed_offset':offset,'arm':arm,'row_libraries_exact':True,'null_DE_calls':int(np.sum(dd['eligible']&(dd['q']<.05))),'actual_DE_calls':int(np.sum(truth['eligible']&(truth['q']<.05))),**diag,**metric,'arrays_sha256':m.b.sha256(arrays)};write(file,rr);del count;gc.collect()
     oo[arm].append((profile,jack,tp,tj,rr));records.append(rr);print('V38',stage,name,i,offset,arm,'fid',rr['fidelity'],flush=True)
   if not null:
    for arm,values in oo.items():
     p,j,tp,tj=(np.array([o[k] for o in values]) for k in range(4));rr=[o[4] for o in values];panel=panel_metrics(p,tp,cc,j,tj,cj,genes,[r['target'] for r in rr],targets)
     for metric in ('fidelity','reach','jaccard','nmae'):
      vv=[r[metric] for r in rr if isinstance(r.get(metric),(float,int)) and np.isfinite(r[metric])];panel[metric]=float(np.mean(vv)) if vv else None
     panel={'family':name,'study':key(name),'seed_offset':offset,'arm':arm,**panel};panels.append(panel);write(dest/f'{name}__{offset}__{arm}__panel.json',panel)
 result={'panels':panels,'records':records,'arms':list(arms),'scope':'Unscaled public diagnostics. No private reference outcomes or anchor scores.'};write(m.A/(stage+'_results.json'),result)
 if not null:
  summary={a:{k:aggregate(panels,a,k) for k in METRICS} for a in arms};write(m.A/(stage+'_summary.json'),summary);print(json.dumps(summary,indent=2),flush=True)
def select():
 check();r=m.read(m.A/'development_results.json');summary=m.read(m.A/'development_summary.json');null=m.read(m.A/'null_results.json')['records'];base=summary['v24'];checks={};scores={}
 for a in r['arms']:
  if a=='v24':continue
  s=summary[a];resolved=all(s.get(k) is not None and base.get(k) is not None for k in METRICS)
  tests={'resolved':resolved,'pds':resolved and s['pds']>=base['pds']-.01,'expression':resolved and s['common_expression_ratio']<=base['common_expression_ratio']*1.02,'fid':resolved and s['fidelity']>=base['fidelity']-.01,'reach':resolved and s['reach']>=base['reach']-.005,'jac':resolved and s['jaccard']>=base['jaccard']-.001,'nmae':resolved and s['nmae']<=base['nmae']*1.01}
  tests['null_no_worse']=all(n['null_DE_calls']<=next(b['null_DE_calls'] for b in null if b['arm']=='v24' and b['family']==n['family']) for n in null if n['arm']==a)
  scores[a]=float(np.mean([(s[k]-base[k])/max(abs(base[k]),.01)*( -1 if k in ('nmae','common_expression_ratio') else 1) for k in METRICS])) if resolved else -100
  tests['positive_average_gain']=scores[a]>0;checks[a]=tests
 eligible=[a for a,tests in checks.items() if all(tests.values())];chosen=max(eligible,key=lambda a:scores[a]) if eligible else 'v24'
 write(m.A/'selection.json',{'selected':chosen,'eligible':eligible,'checks':checks,'relative_gain':scores,'confirmed_improvement':False,'summary':summary});print(json.dumps(m.read(m.A/'selection.json'),indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train','prepare','means','null','development','select','confirmation']);p.add_argument('--outer',choices=list(m.OUTERS)+['all']);p.add_argument('--version',type=int,choices=[38,39],default=38);args=p.parse_args()
 if args.version==39:
  import v39_kolf_hurdle as extra
  extra.install()
 with threadpool_limits(limits=2):
  if args.stage=='train':
   check()
   for outer in [args.outer] if args.outer else m.OUTERS:m.fit(outer)
  elif args.stage in ['null','development','confirmation']:counts(args.stage)
  else:globals()[args.stage]()






