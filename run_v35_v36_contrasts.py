"""Same public cases: two fixed quarter-strength V35 response contrasts."""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v35_methods as m
import v36_methods as pop
from v35_contrast_methods import transform,center_delta
from run_v35_v36 import contexts,key,clean
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
A=m.P36/'contrast';ARMS=['target_centered25','retrieval_contrast25']

def freeze():
 # Mathematical boundary checks before any contrast outcome.
 rng=np.random.default_rng(37);base=rng.uniform(-2.9,2.9,(13,17));pred=np.clip(base+rng.normal(0,.2,base.shape),-3,3);mask=rng.uniform(size=base.shape)>.2;d=center_delta(pred,base,mask);assert abs(d).max()<=.5000001;np.testing.assert_allclose(d.sum(0),0,atol=1e-7);assert (d[~mask]==0).all();assert (base+d<=3+1e-8).all() and (base+d>=-3-1e-8).all()
 A.mkdir(exist_ok=True);assert not (A/'plan.json').exists();m.write(A/'plan.json',{'parent_plan_sha256':m.sha(m.A/'plan.json'),'dependencies':{n:m.sha(m.ROOT/n) for n in ['run_v35_v36_contrasts.py','v35_contrast_methods.py']},'models':{o:m.sha(m.A/'models'/(o+'.npz')) for o in m.OUTERS},'arms':ARMS,'strength':.25,'added_after':'All primary summary/mean/count outcomes and the common-shift diagnostic were known. Explicit exploratory corrective test, not blind validation.','case_rule':'All44 original raw mean cases and identical first4 count cases/context; complete300 predicted targets used for centering, never center only scored cases. Same inputs, seeds and emitter.','hypothesis':'Protect target discrimination: target-centered residual removes target-shared drift; genuine-minus-shuffled residual retains response-retrieval dependence. Fixed quarter strength, no grid.','select_rule':'Prefer a contrast passing expression<=1.02*V24,FID>=V24-.02,PDS>=V24-.02,NMAE<=1.01*V24 and mean MSE/MAE<=1.002, with at least one of mean MSE/MAE/DE metrics improved. Select smallest raw mean MSE among passing contrasts and a currently selected candidate if it passes the same count/mean limits. If none passes, retain current exploratory selection. Strict primary improvement/confirmation flags unchanged.','boundary_checks_passed':True});print('Contrastive corrective test frozen',flush=True)

def run():
 from evaluate_atlas_transfer import load_atlas_model
 m.setup();m.check();plan=json.loads((A/'plan.json').read_text())
 for n,h in plan['dependencies'].items():assert m.sha(m.ROOT/n)==h
 for o,h in plan['models'].items():assert m.sha(m.A/'models'/(o+'.npz'))==h
 assert not (A/'results.json').exists();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[];panels=[];means=[];arrays={}
 for ctx in contexts():
  outer=key(ctx['family']);ex=[outer];x,meta,measured,source=m.raw(ctx['family'],genes);_,cp,bulk,_=moments(x[ctx['left']]);query=cp*1e6;query[~measured]=np.nan;base=m.base_for(query,ex);net,g,rec=m.load(outer,35);pred,_,trace=m.predict(net,g,query,ex,base=base);shuffled,_,_=m.predict(net,g,query,ex,base=base,shuffle=True);effects={a:transform(pred,base,shuffled,genes,targets,g['axis'],trace,query,a) for a in ARMS};ref=x[ctx['right']];refmean=moments(ref)[1]*1e6;outputs={a:[] for a in ARMS}
  for cell in ctx['cells']:
   t=cell['target'];ti=int(np.flatnonzero(targets==t)[0]);own=int(np.flatnonzero(genes==t)[0]);actual=x[cell['actual']];am=moments(actual)[1]*1e6;valid=measured&(refmean>5);valid[own]=False;truth=np.clip(np.log2((am+.1)/(refmean+.1)),-3,3);scores={}
   for arm,eff in {'v24':base,**effects}.items():
    pm=np.nan_to_num(query)*np.exp2(eff[ti]);pm[own]*=.1;pm*=1e6/pm.sum();e=np.clip(np.log2((pm+.1)/(refmean+.1)),-3,3)[valid]-truth[valid];scores[arm]={'square':float(e@e),'absolute':float(abs(e).sum()),'cpm':float(np.sum((pm[valid]-am[valid])**2))}
   for arm in ARMS:means.append({'family':ctx['family'],'study':outer,'target':t,'arm':arm,'mse_ratio':scores[arm]['square']/scores['v24']['square'],'mae_ratio':scores[arm]['absolute']/scores['v24']['absolute'],'cpm_mse_ratio':scores[arm]['cpm']/scores['v24']['cpm']})
  for cell in ctx['cells'][:4]:
   t=cell['target'];ti=int(np.flatnonzero(targets==t)[0]);own=int(np.flatnonzero(genes==t)[0]);template=x[cell['template']];actual=x[cell['actual']];libs=template.sum(1).A1.astype(np.int64);truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual)
   for arm in ARMS:
    count,diag=pop.emit(template,effects[arm][ti],own,cell['seed']+ti,(libs,cp,bulk),atlas);np.testing.assert_array_equal(count.sum(1).A1,libs);dd=de(count,ref,genes);pp,pj=profile_jackknife_chunked(count);outputs[arm].append({'pred':pp,'jack':pj,'truth':tp,'truth_jack':tj});r={'family':ctx['family'],'study':outer,'target':t,'arm':arm,**diag,**direction_metrics(dd,truth,genes,own)};records.append(r)
    for k,v in [('pred_q',dd['q']),('pred_lfc',dd['lfc']),('eligible',truth['eligible']),('truth_q',truth['q']),('truth_lfc',truth['lfc'])]:arrays[ctx['family']+'__'+arm+'__'+t+'__'+k]=v
    print('Contrast',ctx['family'],t,arm,'FID',r['fidelity'],flush=True)
  cc,cj=profile_jackknife_chunked(ref)
  for arm,out in outputs.items():
   pp=np.array([v['pred'] for v in out]);pj=np.array([v['jack'] for v in out]);tp=np.array([v['truth'] for v in out]);tj=np.array([v['truth_jack'] for v in out]);metrics=panel_metrics(pp,tp,cc,pj,tj,cj,genes,[v['target'] for v in ctx['cells'][:4]],targets);rr=[r for r in records if r['family']==ctx['family'] and r['arm']==arm];avg={k:float(np.mean([r[k] for r in rr if np.isfinite(r.get(k,np.nan))])) for k in ['fidelity','reach','jaccard','nmae'] if any(np.isfinite(r.get(k,np.nan)) for r in rr)};panels.append({'family':ctx['family'],'study':outer,'arm':arm,**avg,**metrics})
   for k,v in [('pred',pp),('jack',pj),('truth',tp),('truth_jack',tj),('control',cc),('control_jack',cj)]:arrays[ctx['family']+'__'+arm+'__'+k]=v
  del x;gc.collect()
 summary={a:{k:float(np.mean([np.mean([v[k] for v in means if v['arm']==a and v['study']==study]) for study in m.OUTERS])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']} for a in ARMS};file=A/'audit_arrays.npz';np.savez_compressed(file,**arrays);m.write(A/'results.json',clean({'plan_sha256':m.sha(A/'plan.json'),'records':records,'panels':panels,'mean_records':means,'mean_summary':summary,'audit_arrays_sha256':m.sha(file),'scope':'Exploratory corrective comparison on reused public cases; no fresh confirmation.'}));print(json.dumps(summary,indent=2),flush=True)

def select():
 selection=json.loads((m.P36/'selection.json').read_text());original=m.P36/'precontrast_selection.json';assert not original.exists();original.write_bytes((m.P36/'selection.json').read_bytes());result=json.loads((A/'results.json').read_text());base=selection['count_summaries']['v24'];metrics=['fidelity','reach','jaccard','nmae','pds','common_expression_ratio'];summaries={};checks={};eligible=[]
 for arm in ARMS:
  s={k:float(np.mean([np.mean([r[k] for r in result['panels'] if r['arm']==arm and r['study']==study]) for study in m.OUTERS])) for k in metrics};summaries[arm]=s;mean=result['mean_summary'][arm];guard=s['common_expression_ratio']<=1.02*base['common_expression_ratio'] and s['fidelity']>=base['fidelity']-.02 and s['pds']>=base['pds']-.02 and s['nmae']<=1.01*base['nmae'];improved=mean['mse_ratio']<1 or mean['mae_ratio']<1 or any(s[k]>base[k] for k in ['fidelity','reach','jaccard']);checks[arm]={'count_guard':bool(guard),'mean_noninferiority':mean['mse_ratio']<=1.002 and mean['mae_ratio']<=1.002,'some_gain':bool(improved)}
  if all(checks[arm].values()):eligible.append(arm)
 # Current primary/component candidate enters only if it meets the same limits.
 current=selection['selected_arm'];c=selection['count_summaries'][current];cm=selection['raw_mean_summaries']['v35' if selection['model_version']==35 else 'v36'];current_ok=c['common_expression_ratio']<=1.02*base['common_expression_ratio'] and c['fidelity']>=base['fidelity']-.02 and c['pds']>=base['pds']-.02 and c['nmae']<=1.01*base['nmae'] and cm['mse_ratio']<=1.002 and cm['mae_ratio']<=1.002
 if current_ok:eligible.append(current)
 choice=min(eligible,key=lambda a:result['mean_summary'][a]['mse_ratio'] if a in ARMS else cm['mse_ratio']) if eligible else current
 selection['count_summaries'].update(summaries);selection['contrastive_correction']={'plan_sha256':m.sha(A/'plan.json'),'results_sha256':m.sha(A/'results.json'),'prior_selection_sha256':m.sha(original),'checks':checks,'eligible':eligible,'selected':choice,'exploratory_after_primary_outcomes':True};selection['raw_mean_summaries'].update(result['mean_summary'])
 if choice in ARMS:selection['selected_arm']=choice;selection['model_version']=35;selection.pop('dispersion_model_version',None)
 m.write(m.P36/'selection.json',selection);m.write(A/'decision.json',selection['contrastive_correction']);print(json.dumps({'selected':choice,'checks':checks,'means':result['mean_summary'],'counts':summaries},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run','select']);args=p.parse_args()
 with threadpool_limits(limits=2):{'freeze':freeze,'run':run,'select':select}[args.stage]()