"""Exploratory component check: V35 means with V36 dispersion only.

Added after mean/partial count outcomes were visible; same cases and fixed
weights. Never presented as a new blind or preregistered primary comparison.
"""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v35_methods as m
import v36_methods as p
from run_v35_v36 import contexts,key,clean
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
A=m.P36/'hybrid';ARM='v35_population_dispersion'

def freeze():
 A.mkdir(exist_ok=True);assert not (A/'plan.json').exists()
 models={str(f.relative_to(m.ROOT)):m.sha(f) for folder in [m.A,m.P36] for f in (folder/'models').glob('*.npz') if f.stem!='all'}
 m.write(A/'plan.json',{'code_sha256':m.sha(m.ROOT/'run_v35_v36_hybrid.py'),'parent_plan_sha256':m.sha(m.A/'plan.json'),'models':models,'arm':ARM,'added_after':'V35/V36 summary and44 raw mean results, plus partial primary treatment-count outcomes. Explicit exploratory extension, not new confirmation.','hypothesis':'Keep stronger V35 means; apply only the fixed V36 treatment-specific variance head. Head is conditioned on V24/source evidence, never on an observed treatment or predicted V36 mean.','cases':'Identical first4 cases/context, same templates, controls, source exclusions and seeds. No new blend or hyperparameter search.','adopt_if':'Only when the primary selector chooses V35: hybrid passes expression<=1.02*V24,FID>=V24-.02,PDS>=V24-.02; improves>=2 of FID/REACH/JAC/NMAE versus V35, and expression error does not increase versus V35. Null behavior is identical by construction, as dispersion/effect zero without a perturbation. All strict primary improvement/confirmation flags retained.'});print('Exploratory hybrid plan frozen',flush=True)

def run():
 from evaluate_atlas_transfer import load_atlas_model
 m.setup();m.check();plan=json.loads((A/'plan.json').read_text());assert plan['code_sha256']==m.sha(m.ROOT/'run_v35_v36_hybrid.py')
 for f,h in plan['models'].items():assert m.sha(m.ROOT/f)==h
 assert not (A/'results.json').exists();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[];panels=[];arrays={}
 for ctx in contexts():
  outer=key(ctx['family']);ex=[outer];x,meta,measured,source=m.raw(ctx['family'],genes);ctrl=m.control_stats(x[ctx['left']]);_,cp,bulk,_=moments(x[ctx['left']]);query=cp*1e6;query[~measured]=np.nan;base=m.base_for(query,ex);cells=ctx['cells'][:4];ts=np.array([c['target'] for c in cells]);idx=[int(np.flatnonzero(targets==t)[0]) for t in ts];net,g,_=m.load(outer,35);effects,_,_=m.predict(net,g,query,ex,ts,base[idx]);net,g,_=m.load(outer,36);_,disp,_=m.predict(net,g,query,ex,ts,base[idx],control=ctrl);ref=x[ctx['right']];out=[]
  for i,cell in enumerate(cells):
   t=cell['target'];own=int(np.flatnonzero(genes==t)[0]);template=x[cell['template']];actual=x[cell['actual']];libs=template.sum(1).A1.astype(np.int64);count,diag=p.emit(template,effects[i],own,cell['seed']+idx[i],(libs,cp,bulk),atlas,disp[i]);np.testing.assert_array_equal(count.sum(1).A1,libs);truth=de(actual,ref,genes);pred=de(count,ref,genes);pp,pj=profile_jackknife_chunked(count);tp,tj=profile_jackknife_chunked(actual);out.append({'pred':pp,'jack':pj,'truth':tp,'truth_jack':tj});r={'family':ctx['family'],'study':outer,'target':t,'arm':ARM,**diag,**direction_metrics(pred,truth,genes,own)};records.append(r)
   for name,value in [('pred_q',pred['q']),('pred_lfc',pred['lfc']),('eligible',truth['eligible']),('truth_q',truth['q']),('truth_lfc',truth['lfc'])]:arrays[ctx['family']+'__'+t+'__'+name]=value
   print('Hybrid',ctx['family'],t,'FID',r['fidelity'],flush=True)
  cc,cj=profile_jackknife_chunked(ref);pp=np.array([v['pred'] for v in out]);pj=np.array([v['jack'] for v in out]);tp=np.array([v['truth'] for v in out]);tj=np.array([v['truth_jack'] for v in out]);metrics=panel_metrics(pp,tp,cc,pj,tj,cj,genes,ts,targets);rr=[r for r in records if r['family']==ctx['family']];avg={k:float(np.mean([r[k] for r in rr if np.isfinite(r.get(k,np.nan))])) for k in ['fidelity','reach','jaccard','nmae'] if any(np.isfinite(r.get(k,np.nan)) for r in rr)};panels.append({'family':ctx['family'],'study':outer,'arm':ARM,**avg,**metrics})
  for name,value in [('pred',pp),('jack',pj),('truth',tp),('truth_jack',tj),('control',cc),('control_jack',cj)]:arrays[ctx['family']+'__'+name]=value
  del x;gc.collect()
 file=A/'audit_arrays.npz';np.savez_compressed(file,**arrays);m.write(A/'results.json',clean({'plan_sha256':m.sha(A/'plan.json'),'records':records,'panels':panels,'audit_arrays_sha256':m.sha(file),'scope':'Exploratory component test on the same spent20 public cases. No new blind evidence.'}))

def select():
 selection=json.loads((m.P36/'selection.json').read_text());assert 'hybrid_component' not in selection;original=m.P36/'primary_selection.json';assert not original.exists();original.write_bytes((m.P36/'selection.json').read_bytes());result=json.loads((A/'results.json').read_text());metrics=['fidelity','reach','jaccard','nmae','pds','common_expression_ratio'];summary={}
 for k in metrics:
  vals=[]
  for study in m.OUTERS:
   rr=[r[k] for r in result['panels'] if r['study']==study and isinstance(r.get(k),(float,int))];vals.append(float(np.mean(rr)) if rr else np.nan)
  summary[k]=float(np.mean(vals)) if np.isfinite(vals).all() else None
 base=selection['count_summaries']['v24'];ref=selection['count_summaries']['v35'];finite=all(summary[k] is not None for k in metrics);guard=finite and summary['common_expression_ratio']<=1.02*base['common_expression_ratio'] and summary['fidelity']>=base['fidelity']-.02 and summary['pds']>=base['pds']-.02;better=sum(summary[k] is not None and ref[k] is not None and (summary[k]<ref[k] if k=='nmae' else summary[k]>ref[k]) for k in ['fidelity','reach','jaccard','nmae']);adopt=selection['selected_arm']=='v35' and guard and better>=2 and summary['common_expression_ratio']<=ref['common_expression_ratio'];selection['count_summaries'][ARM]=summary;selection['hybrid_component']={'plan_sha256':m.sha(A/'plan.json'),'results_sha256':m.sha(A/'results.json'),'primary_selection_sha256':m.sha(original),'count_guard_passes':bool(guard),'DE_metrics_improved_vs_v35':int(better),'adopted':bool(adopt),'exploratory_after_partial_outcomes':True}
 if adopt:selection['selected_arm']=ARM;selection['dispersion_model_version']=36
 m.write(m.P36/'selection.json',selection);m.write(A/'decision.json',selection['hybrid_component']);print(json.dumps({'selected_arm':selection['selected_arm'],'hybrid':selection['hybrid_component'],'summary':summary},indent=2),flush=True)
if __name__=='__main__':
 p0=argparse.ArgumentParser();p0.add_argument('stage',choices=['freeze','run','select']);args=p0.parse_args()
 with threadpool_limits(limits=2):{'freeze':freeze,'run':run,'select':select}[args.stage]()