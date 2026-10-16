"""Sequential V35 retrieval, then V36 population-view experiments and raw counts."""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v35_methods as m
import v36_methods as pop
from run_v30_experiment import rows_for,summarize_context
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
from v6_methods import moments


def clean(v):
 if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
 if isinstance(v,list):return [clean(x) for x in v]
 if isinstance(v,(np.integer,)):return int(v)
 if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else None
 return v

def query_moments(source):
 if source not in m.RAW_NAME:return None
 d=m.old.old.read_cache();x,meta,measured,rec=m.raw(m.RAW_NAME[source],d['genes']);ids=np.flatnonzero(meta.gene_target.to_numpy(str)=='Non-Targeting');v=m.control_stats(x[ids]);del x;gc.collect();return v

def summary(version):
 m.setup();m.check();folder=m.A if version==35 else m.P36;assert not (folder/'summary_results.json').exists();data=m.old.old.read_cache();records=[];predfolder=folder/'predictions';predfolder.mkdir(exist_ok=True)
 for outer in m.OUTERS:
  net,g,model=m.load(outer,version);rr=rows_for(data,outer);ex=[outer]
  for ci in np.unique(data['row_contexts'][rr]):
   rows=rr[data['row_contexts'][rr]==ci];query=data['profiles'][ci];source=next(s for s,v in m.CI.items() if v==ci);control=query_moments(source) if version==36 else None;base=m.base_for(query,ex);save={'base':base,'query':query,'rows':rows}
   for arm,shuffle,static in [('learned',False,False),('response_shuffled',True,False),('static_retrieval',False,True)]:
    pred,var,trace=m.predict(net,g,query,ex,base=base,shuffle=shuffle,static=static,control=control);out=summarize_context(data,rows,pred,base);records.append({'outer':outer,'context':int(ci),'source':source,'arm':arm,**out});save[arm]=pred;save[arm+'_dispersion']=var
    print('V'+str(version),'summary',outer,source,arm,out['mse_ratio'],out['mae_ratio'],flush=True)
   path=predfolder/(source+'.npz');np.savez_compressed(path,**save);m.write(path.with_suffix('.json'),{'sha256':m.sha(path),'model_sha256':model['sha256'],'outer':outer,'source':source,'excluded':ex,'plan_sha256':m.sha(m.A/'plan.json'),'axis':g['axis'].tolist(),'retrieval_trace':trace})
  del net;gc.collect()
 study={arm:{o:{k:float(np.mean([r[k] for r in records if r['outer']==o and r['arm']==arm])) for k in ['mse_ratio','mae_ratio']} for o in m.OUTERS} for arm in ['learned','response_shuffled','static_retrieval']}
 avg={arm:{k:float(np.mean([v[k] for v in study[arm].values()])) for k in ['mse_ratio','mae_ratio']} for arm in study};s=avg['learned'];checks={'MSE_gain':s['mse_ratio']<=.99,'MAE_gain':s['mae_ratio']<=.995,'worst_family':max(v['mse_ratio'] for v in study['learned'].values())<=1.02,'beats_shuffle':all(s[k]<avg['response_shuffled'][k] for k in s)}
 m.write(folder/'summary_results.json',{'plan_sha256':m.sha(m.A/'plan.json'),'records':records,'study_summary':study,'summary':avg,'checks':checks,'passes':all(checks.values()),'scope':'Reused public summary LFC development, entire recipient family excluded. Corruption/static controls are inference ablations, not separately trained shuffled models.'})
 print(json.dumps({'version':version,'summary':avg,'checks':checks},indent=2),flush=True)

def contexts():return json.loads((m.ROOT/'experiments/v34_residual_flow/analysis/plan.json').read_text())['contexts']
def key(ctx):return 'CD4' if ctx.startswith('CD4') else ('Orion' if ctx in ['HCT116','HEK293T'] else ctx)

def mean_screen():
 m.setup();m.check();records=[];data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];folder=m.P36;assert not (folder/'mean_results.json').exists();models={}
 for ctx in contexts():
  outer=key(ctx['family']);ex=[outer];x,meta,measured,source=m.raw(ctx['family'],genes);controls=x[ctx['left']];ctrl=m.control_stats(controls);query=ctrl[0]*1e6;query[~measured]=np.nan;ref=moments(x[ctx['right']])[1]*1e6;base=m.base_for(query,ex);tlist=np.array([c['target'] for c in ctx['cells']]);indices=np.array([int(np.flatnonzero(targets==t)[0]) for t in tlist]);pred={'v24':base[indices]}
  for ver in [35,36]:
   if (outer,ver) not in models:models[outer,ver]=m.load(outer,ver)
   net,g,_=models[outer,ver];p,v,_=m.predict(net,g,query,ex,tlist,base[indices],control=ctrl if ver==36 else None);pred['v'+str(ver)]=p
  for i,cell in enumerate(ctx['cells']):
   own=int(np.flatnonzero(genes==cell['target'])[0]);actual=moments(x[cell['actual']])[1]*1e6;valid=measured&(ref>5);valid[own]=False;truth=np.clip(np.log2((actual+.1)/(ref+.1)),-3,3);scores={}
   for arm,eff in pred.items():
    p=np.nan_to_num(query)*np.exp2(eff[i]);p[own]*=.1;p*=1e6/p.sum();lfc=np.clip(np.log2((p+.1)/(ref+.1)),-3,3);delta=lfc[valid]-truth[valid];scores[arm]={'square':float(delta@delta),'absolute':float(abs(delta).sum()),'cpm':float(np.sum((p[valid]-actual[valid])**2))}
   for arm,s in scores.items():records.append({'family':ctx['family'],'study':outer,'target':cell['target'],'arm':arm,'mse_ratio':s['square']/scores['v24']['square'],'mae_ratio':s['absolute']/scores['v24']['absolute'],'cpm_mse_ratio':s['cpm']/scores['v24']['cpm']})
  print('V35/V36 raw mean',ctx['family'],len(ctx['cells']),flush=True);del x;gc.collect()
 study={arm:{o:{k:float(np.mean([r[k] for r in records if r['study']==o and r['arm']==arm])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']} for o in m.OUTERS} for arm in ['v24','v35','v36']};summary={a:{k:float(np.mean([v[k] for v in study[a].values()])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']} for a in study}
 m.write(folder/'mean_results.json',{'plan_sha256':m.sha(m.A/'plan.json'),'records':records,'study_summary':study,'summary':summary,'scope':'Same44 previously spent public raw mean cases; no blind confirmation or private anchored score.'});print(json.dumps(summary,indent=2),flush=True)

def count_test(stage):
 from evaluate_atlas_transfer import load_atlas_model
 from run_v26_experiment import distribution
 m.setup();m.check();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[];panels=[];save={};models={};assert not (m.P36/(stage+'_results.json')).exists()
 for ctx in contexts():
  outer=key(ctx['family']);ex=[outer];null=stage=='null';x,meta,measured,source=m.raw(ctx['family'],genes);ids=ctx['null']['fit_controls'] if null else ctx['left'];ref=x[ctx['null']['reference'] if null else ctx['right']];ctrl=m.control_stats(x[ids]);depth,p,bulk,_=moments(x[ids]);query=p*1e6;query[~measured]=np.nan;cells=[ctx['null']] if null else ctx['cells'][:4];out=[]
  if not null:
   base=m.base_for(query,ex);tlist=np.array([c['target'] for c in cells]);indices=np.array([int(np.flatnonzero(targets==t)[0]) for t in tlist]);effects={'v24':base[indices]};dispersions={}
   for ver in [35,36]:
    if (outer,ver) not in models:models[outer,ver]=m.load(outer,ver)
    net,g,_=models[outer,ver];pred,dd,trace=m.predict(net,g,query,ex,tlist,base[indices],control=ctrl if ver==36 else None);effects['v'+str(ver)]=pred;dispersions['v'+str(ver)]=dd
  for ii,cell in enumerate(cells):
   t=cell['target'];own=None if null else int(np.flatnonzero(genes==t)[0]);ti=0 if null else int(np.flatnonzero(targets==t)[0]);seed=cell['seed']+ti;template=x[cell['template']];actual=x[cell['actual']];libs=template.sum(1).A1.astype(np.int64);prepared=(libs,p,bulk);arms={};diagnostic={};zero=np.zeros(len(genes))
   for arm in ['v24','v35','v36','v36_dispersion']:
    eff=zero if null else effects['v36' if arm=='v36_dispersion' else arm][ii];dd=(zero if null else dispersions['v36'][ii]) if arm=='v36_dispersion' else None
    count,diag=pop.emit(template,eff,own,seed,prepared,atlas,dd);arms[arm]=count;diagnostic[arm]=diag
   truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);entry={'truth':tp,'truth_jack':tj};axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));save[ctx['family']+'__'+t+'__truth_q']=truth['q'];save[ctx['family']+'__'+t+'__truth_lfc']=truth['lfc'];save[ctx['family']+'__'+t+'__truth_eligible']=truth['eligible']
   for arm,count in arms.items():
    np.testing.assert_array_equal(count.sum(1).A1,libs);assert count.dtype.kind in 'iu' and (count.data>=0).all();pp,pj=profile_jackknife_chunked(count);entry[arm+'_pred']=pp;entry[arm+'_jack']=pj;d=de(count,ref,genes);metric=direction_metrics(d,truth,genes,-1 if own is None else own);rec={'family':ctx['family'],'study':outer,'target':t,'arm':arm,'predicted_cells':count.shape[0],'actual_cells':actual.shape[0],'null_DE_calls':int(((d['q']<.05)&d['eligible']).sum()),'real_null_DE_calls':int(((truth['q']<.05)&truth['eligible']).sum()),**diagnostic[arm],**metric,**distribution(count,ref,axis)};records.append(rec)
    for k0 in ['q','lfc','eligible']:save[ctx['family']+'__'+t+'__'+arm+'_'+k0]=d[k0]
    print('V35/V36',stage,ctx['family'],t,arm,'FID',rec['fidelity'],'calls',rec['null_DE_calls'],flush=True)
   out.append(entry)
  cp,cj=profile_jackknife_chunked(ref)
  for arm in arms:
   pp=np.array([v[arm+'_pred'] for v in out]);pj=np.array([v[arm+'_jack'] for v in out]);tp=np.array([v['truth'] for v in out]);tj=np.array([v['truth_jack'] for v in out]);metrics={} if null else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[v['target'] for v in cells],targets);rr=[v for v in records if v['family']==ctx['family'] and v['arm']==arm];av={k:float(np.mean([r[k] for r in rr if np.isfinite(r.get(k,np.nan))])) for k in ['fidelity','reach','jaccard','nmae'] if any(np.isfinite(r.get(k,np.nan)) for r in rr)};panels.append({'family':ctx['family'],'study':outer,'arm':arm,**av,**metrics})
   for k0,value in [('pred',pp),('jack',pj),('truth',tp),('truth_jack',tj),('control',cp),('control_jack',cj)]:save[ctx['family']+'__'+arm+'__'+k0]=value
  del x;gc.collect()
 file=m.P36/(stage+'_audit_arrays.npz');np.savez_compressed(file,**save);m.write(m.P36/(stage+'_results.json'),clean({'plan_sha256':m.sha(m.A/'plan.json'),'records':records,'panels':panels,'audit_arrays_sha256':m.sha(file),'scope':'Raw public count diagnostics; fixed20 treatment cases or5 nulls; D1 treatment200cells, remaining400. No official anchor scaling.'}))

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train35','summary35','train36','summary36','mean','null','treatment']);p.add_argument('--outer',choices=m.OUTERS+['all']);a=p.parse_args()
 if a.stage=='freeze':m.freeze()
 elif a.stage in ['train35','train36']:
  for o in [a.outer] if a.outer else m.OUTERS:
   (m.train35 if a.stage=='train35' else pop.train)(o)
 elif a.stage.startswith('summary'):summary(int(a.stage[-2:]))
 elif a.stage=='mean':mean_screen()
 else:count_test(a.stage)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()