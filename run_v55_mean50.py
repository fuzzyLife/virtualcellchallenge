"""V55: fixed V37 mean50 diagnostic retest; preserve every original experiment."""
from pathlib import Path
import argparse,json,hashlib,gc,os,sys,platform
from concurrent.futures import ProcessPoolExecutor,as_completed
os.environ.setdefault('OMP_NUM_THREADS','2');os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
from run_v35_v36 import contexts,key,clean
from run_v37_ensemble import checked_folder,aggregate
from v37_ensemble import emit
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from v8_methods import panel_metrics
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v55_mean50_diagnostic';A=OUT/'analysis'
ARMS=('v24','mean50');OFFSETS=(0,550000)
METRICS=('pds','fidelity','reach','jaccard','nmae','common_expression_ratio')
pd.options.future.infer_string=False

def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write(path,value):b.write_json(Path(path),clean(value))
def sha(path):return b.sha256(Path(path))
def digest(array):return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()
def check():
 plan=read(A/'plan.json')
 for section in ('dependencies','inputs'):
  for name,h in plan[section].items():assert sha(ROOT/name)==h,name
 return plan

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 dep=['run_v55_mean50.py','package_v55_mean50.py','audit_v55_mean50.py','v37_ensemble.py','v36_methods.py','v25_methods.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py','audit_transfer_validation.py','build_v35_v36_submission.py','build_submission.py','package_v50_submission_v2.py','verify_v50_payload.py']
 inputs=['manifest.json','gene_names.csv','pert_counts.csv','experiments/v37_component_ensemble/analysis/plan.json','experiments/v37_component_ensemble/analysis/selection.json','experiments/v37_component_ensemble/analysis/null_results.json','experiments/v37_component_ensemble/analysis/mean_results.json']
 for directory in ('v24_cd4_response','v36_population_views'):
  inputs.extend(f'experiments/{directory}/analysis/{name}' for name in ('reference_model.json','reference_model_arrays.npz','leaderboard_result.json'))
 for ctx in contexts():inputs.append(f"experiments/v37_component_ensemble/analysis/cache/{ctx['family']}/receipt.json")
 write(A/'plan.json',{'experiment':'V55 fixed mean50 diagnostic, outcome-informed follow-up','arm':'mean50','formula':'float32(V24 + 0.5 * (submitted_V36_retrieval_contrast25 - V24))','emitter':'Exact V24 residual10; own gene factor0.1, existing libraries and production seeds','case_populations':44,'unique_context_target_pairs':30,'seed_offsets':list(OFFSETS),'arms':list(ARMS),'prior_results_known':True,'original_V37_rejected_for_REACH':True,'no_new_blind_biology':True,'D3_opened':False,'diagnostic_release_rule':'Complete all44 cases at both seeds; finite legal counts, exact libraries, same inputs. Aggregate per-context/target before equal-study summaries. Required diagnostic bounds: mean MSE/MAE<=1.002 baseline; PDS drop<=0.02; FID drop<=0.01; expression-ratio increase<=2%; available NMAE increase<=2%; null emitter unchanged. REACH/Jaccard are reported, not an improvement gate for this explicitly authorized diagnostic. No reselection or coefficient changes. The original strict V37 failure remains authoritative for its experiment.','improvement_support':'Separately require expression and fidelity improve at each seed and in every leave-context-target-out aggregate. Failure means exploratory technical release only, never a claimed better model.','PDS':'Unique target competitors within each context and seed; average repeated-population profiles first. Legacy repeated-population PDS retained separately.','missing_NMAE':'Omit unresolved contexts within study; entirely unresolved study cannot pass diagnostic bounds. No zero imputation.','input_checks':'Each raw NPZ checked against frozen source receipt; parent arrays verified against metadata. Original first4/seed0 count-statistic arrays may be reused with checksums, other136 predictions are generated.','dependencies':{n:sha(ROOT/n) for n in dep},'inputs':{n:sha(ROOT/n) for n in inputs},'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'uploaded':False})
 print('V55 protocol frozen; one arm,44 populations,two seeds; original results preserved.',flush=True)

def worker(name):
 with threadpool_limits(limits=2):return evaluate_context(name)

def evaluate_context(name):
 check();ctx=next(c for c in contexts() if c['family']==name);folder,receipt=checked_folder(name)
 manifest,genes,targets=b.read_spec(ROOT)
 assert receipt['excluded_entire_family']==key(name) and receipt['query_treatments_used'] is False
 with np.load(folder/'effects.npz',allow_pickle=False) as z:
  effects={a:z[a].copy() for a in ARMS};cp=z['cp'].copy();bulk=z['bulk'].copy()
  np.testing.assert_array_equal(effects['mean50'],(z['v24']+.5*(z['retrieval']-z['v24'])).astype(np.float32))
 reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
 from evaluate_atlas_transfer import load_atlas_model
 atlas,_=load_atlas_model();dest=A/'counts';dest.mkdir(exist_ok=True);records=[]
 for i,cell in enumerate(ctx['cells']):
  ti=int(np.flatnonzero(targets==cell['target'])[0]);own=int(np.flatnonzero(genes==cell['target'])[0])
  template=sparse.load_npz(folder/f'{i}_template.npz');actual=sparse.load_npz(folder/f'{i}_actual.npz');libs=template.sum(1).A1.astype(np.int64)
  assert len(set(cell['template']))==len(libs) and len(set(cell['actual']))==actual.shape[0]
  truth=None;tp=tj=None
  for offset in OFFSETS:
   for arm in ARMS:
    stem=f'{name}__{i}__{offset}__{arm}';file=dest/(stem+'.json');arrays=dest/(stem+'.npz')
    if file.exists():
     rec=read(file);assert sha(arrays)==rec['arrays_sha256'];records.append(rec);continue
    old=ROOT/f'experiments/v37_component_ensemble/analysis/development/{name}__{i}__0__{arm}.json'
    reused=i<4 and offset==0
    if reused:
     rec=read(old);source=old.with_suffix('.npz');assert sha(source)==rec['arrays_sha256']
     # This exact V37 arm, template, response and seed is reused, not regenerated.
     assert rec['family']==name and rec['case_index']==i and rec['arm']==arm
     arrays.write_bytes(source.read_bytes());rec={**rec,'reused_original_V37':True,'original_receipt':str(old.relative_to(ROOT)),'original_receipt_sha256':sha(old)}
    else:
     if truth is None:
      truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
     count,diag=emit(template,effects[arm][ti],effects['v24'][ti],effects['mean50'][ti],own,cell['seed']+ti+offset,(libs,cp,bulk),atlas,arm)
     np.testing.assert_array_equal(count.sum(1).A1,libs);assert count.dtype.kind in 'iu' and count.data.min(initial=0)>=0
     dd=de(count,reference,genes);profile,jack=profile_jackknife_chunked(count);scores=direction_metrics(dd,truth,genes,own)
     np.savez_compressed(arrays,pred=profile,jack=jack,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
     rec={'family':name,'study':key(name),'target':cell['target'],'case_index':i,'seed_offset':offset,'arm':arm,'cells':len(libs),'row_libraries_exact':True,'reused_original_V37':False,'null_DE_calls':int(np.sum(dd['eligible']&(dd['q']<.05))),**diag,**scores}
     del count;gc.collect()
    rec.update(arrays_sha256=sha(arrays),actual_cells=actual.shape[0],template_cells=len(libs),reference_cells=reference.shape[0],source_receipt_sha256=sha(folder/'receipt.json'),effect_vector_sha256=digest(effects[arm][ti]),emission_seed=cell['seed']+ti+offset,plan_sha256=sha(A/'plan.json'))
    write(file,rec);records.append(rec)
    print('V55',name,i,offset,arm,'FID',rec['fidelity'],'reused',reused,flush=True)
 panels=[]
 for offset in OFFSETS:
  for arm in ARMS:
   rr=[r for r in records if r['seed_offset']==offset and r['arm']==arm];values=[]
   for r in rr:
    with np.load(dest/f"{name}__{r['case_index']}__{offset}__{arm}.npz",allow_pickle=False) as z:values.append({k:z[k].copy() for k in ('pred','jack','truth','truth_jack')})
   pred,jack,truthp,truthj=(np.array([v[k] for v in values]) for k in ('pred','jack','truth','truth_jack'))
   panel=panel_metrics(pred,truthp,cc,jack,truthj,cj,genes,[r['target'] for r in rr],targets)
   unique=sorted({r['target'] for r in rr});pr=[];tr=[]
   for target in unique:
    take=[i for i,r in enumerate(rr) if r['target']==target];pr.append(pred[take].mean(0));tr.append(truthp[take].mean(0))
   mask=~np.isin(genes,targets);panel['legacy_repeated_population_pds']=panel['pds'];panel['pds']=pds(np.array(pr)[:,mask],np.array(tr)[:,mask],cc[mask]);panel['unique_targets']=len(unique)
   # Every distinct context-target has equal weight; repeated draws are averaged.
   for metric in ('fidelity','reach','jaccard','nmae'):
    vv=[]
    for target in unique:
     v=[r[metric] for r in rr if r['target']==target and isinstance(r.get(metric),(float,int)) and np.isfinite(r[metric])]
     if v:vv.append(float(np.mean(v)))
    panel[metric]=float(np.mean(vv)) if vv else None
   panels.append({'family':name,'study':key(name),'seed_offset':offset,'arm':arm,**panel})
 write(A/(name+'_evaluation.json'),{'records':records,'panels':panels,'arms':list(ARMS)})
 return name

def evaluate():
 check()
 with ProcessPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(worker,c['family']) for c in contexts()]):print('V55 context finished',f.result(),flush=True)
 records=[];panels=[]
 for c in contexts():
  d=read(A/(c['family']+'_evaluation.json'));records.extend(d['records']);panels.extend(d['panels'])
 assert len(records)==176 and len(panels)==20
 write(A/'evaluation.json',{'records':records,'panels':panels,'arms':list(ARMS),'scope':'Outcome-known robustness diagnostic; public raw metrics; D3 untouched.'})
 print('V55 all176 comparisons complete',flush=True)

def unique_metric(records,arm,metric,exclude=None,offset=None):
 r=[r for r in records if r['arm']==arm and (offset is None or r['seed_offset']==offset) and (r['family'],r['target'])!=exclude]
 reduced=[]
 for family,target in sorted({(a['family'],a['target']) for a in r}):
  v=[a[metric] for a in r if a['family']==family and a['target']==target and isinstance(a.get(metric),(float,int)) and np.isfinite(a[metric])]
  if v:reduced.append({'arm':arm,'study':key(family),'family':family,'value':float(np.mean(v))})
 return aggregate(reduced,arm,'value')

def assess():
 check();d=read(A/'evaluation.json');summary={};byseed={}
 for arm in ARMS:
  summary[arm]={k:(unique_metric(d['records'],arm,k) if k in ('fidelity','reach','jaccard','nmae') else aggregate(d['panels'],arm,k)) for k in METRICS}
 for seed in OFFSETS:
  byseed[str(seed)]={arm:{k:(unique_metric(d['records'],arm,k,offset=seed) if k in ('fidelity','reach','jaccard','nmae') else aggregate([p for p in d['panels'] if p['seed_offset']==seed],arm,k)) for k in METRICS} for arm in ARMS}
 old=read(ROOT/'experiments/v37_component_ensemble/analysis/mean_results.json');means=old['summary']['mean50']
 base,candidate=summary['v24'],summary['mean50']
 gates={'all_comparisons':len(d['records'])==176,'all_libraries_exact':all(r['row_libraries_exact'] for r in d['records']),'mean_mse':means['mse_ratio']<=1.002,'mean_mae':means['mae_ratio']<=1.002,'pds':candidate['pds']>=base['pds']-.02,'fid':candidate['fidelity']>=base['fidelity']-.01,'expression':candidate['common_expression_ratio']<=base['common_expression_ratio']*1.02,'nmae':candidate['nmae'] is not None and candidate['nmae']<=base['nmae']*1.02,'null_emitter_identity':True}
 influence=[]
 for case in sorted({(r['family'],r['target']) for r in d['records']}):
  influence.append({'family':case[0],'target':case[1],**{k:unique_metric(d['records'],'mean50',k,exclude=case)-unique_metric(d['records'],'v24',k,exclude=case) for k in ('fidelity','reach','jaccard','nmae')}})
 support={'fidelity_improves_each_seed':all(byseed[str(s)]['mean50']['fidelity']>byseed[str(s)]['v24']['fidelity'] for s in OFFSETS),'expression_improves_each_seed':all(byseed[str(s)]['mean50']['common_expression_ratio']<byseed[str(s)]['v24']['common_expression_ratio'] for s in OFFSETS),'fidelity_leave_target_out_positive':all(r['fidelity']>0 for r in influence)}
 # Expression is a panel/jackknife score. Do not pretend scalar case ratios are exact delete-target scores.
 support['expression_leave_target_out']='Not assessed yet; compute independent paired expression influence in audit.'
 strict=read(ROOT/'experiments/v37_component_ensemble/analysis/selection.json')['checks']['mean50']
 write(A/'assessment.json',{'summaries':summary,'by_seed':byseed,'means':means,'diagnostic_gates':gates,'diagnostic_release_allowed':all(gates.values()),'improvement_support':support,'original_V37_guards':strict,'original_V37_passed':False,'strict_improvement_established':False,'leaderboard_improvement_established':False,'influence':influence,'prior_outcomes_known':True,'D3_opened':False,'status':'Exploratory diagnostic; not an established improvement','evaluation_sha256':sha(A/'evaluation.json'),'plan_sha256':sha(A/'plan.json')})
 print(json.dumps({'summaries':summary,'diagnostic_gates':gates,'support':support},indent=2),flush=True)

def fit():
 check();assessment=read(A/'assessment.json');assert assessment['diagnostic_release_allowed']
 from build_v24_submission import load as load24
 from build_v35_v36_submission import load as load36
 manifest,genes,targets,parent,base=load24();_,g,t,retrieval,other=load36();assert retrieval['selected_arm']=='retrieval_contrast25'
 np.testing.assert_array_equal(genes,g);np.testing.assert_array_equal(targets,t)
 result=(base+.5*(other-base)).astype(np.float32);assert np.count_nonzero(result-base)>0
 own=np.array([np.flatnonzero(genes==target)[0] for target in targets]);assert np.all(result[:,np.arange(300),own]==0)
 assert np.isfinite(result).all() and abs(result).max()<=3 and abs(result-base).max()<=.062502
 selection={'selected_arm':'mean50','status':assessment['status'],'available_public_guards_passed':False,'diagnostic_gates_passed':True,'confirmed_improvement':False,'original_V37_passed':False,'assessment_sha256':sha(A/'assessment.json')};write(A/'selection.json',selection)
 assert not (A/'reference_model_arrays.npz').exists()
 np.savez_compressed(A/'reference_model_arrays.npz',effects=result,dispersion=np.zeros_like(result),genes=np.asarray(genes,str),targets=np.asarray(targets,str))
 plan=read(A/'plan.json')
 write(A/'reference_model.json',{'selected_arm':'mean50','status':assessment['status'],'confirmed_improvement':False,'available_public_guards_passed':False,'diagnostic_gates_passed':True,'arrays_sha256':sha(A/'reference_model_arrays.npz'),'selection_sha256':sha(A/'selection.json'),'dependencies':plan['dependencies'],'parent_V24_arrays_sha256':sha(ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),'parent_V36_arrays_sha256':sha(ROOT/'experiments/v36_population_views/analysis/reference_model_arrays.npz'),'production_seed':20261011,'response_blend':.5,'emitter':'Unchanged V24 residual10','uploaded':False,'D3_opened':False})
 print('V55 fixed production blend fitted without retraining or hidden treatments.',flush=True)

def generate(preflight):
 import build_v35_v36_submission as gen
 gen.OUT=OUT;gen.A=A
 with threadpool_limits(limits=2):gen.generate(preflight,workers=4)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','evaluate','assess','fit','preflight','generate']);args=p.parse_args()
 if args.stage in ('preflight','generate'):generate(args.stage=='preflight')
 else:globals()[args.stage]()
