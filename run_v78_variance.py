"""Paired V78 shape trials with exact profile invariants and all DE evidence."""
import argparse,gc,json,time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v78_margin_variance as m
from v57_pie import read,write,sha
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
from audit_v76_counts import spot_ranks

def cases(ctx):
 seen=set();out=[]
 for i,c in enumerate(ctx['cells']):
  if c['target'] not in seen:seen.add(c['target']);out.append(i)
 return out

def counts(name,stage,offset):
 m.check();d=m.data();genes=d['genes'];targets=d['targets'];ctx=next(c for c in m.contexts() if c['family']==name);folder,_=m.checked_folder(name)
 if stage=='development':arms=('v66',)+m.ARMS+tuple(a+'_shuffled' for a in m.ARMS);indices=cases(ctx)[:2]
 else:selected=read(m.A/'selection.json')['selected'];assert selected;arms=('v66',selected,selected+'_shuffled');indices=cases(ctx)[2:]
 query=m.query(name);np.testing.assert_array_equal(query['genes'],genes);panel=np.isin(genes,targets)
 genuine,trace=m.bridge.pooled((m.key(name),));fake,faketrace=m.bridge.pooled((m.key(name),),True)
 parent=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==read(parent.with_suffix('.json'))['sha256']
 with np.load(parent,allow_pickle=False) as z:effects=z['assay75'].copy()
 with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy();bulk=z['bulk'].copy()
 atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
 dest=m.A/stage/str(offset);dest.mkdir(parents=True,exist_ok=True);outputs={a:[] for a in arms}
 for i in indices:
  cell=ctx['cells'][i];target=cell['target'];ti=list(targets).index(target);own=list(genes).index(target)
  template=sparse.load_npz(folder/f'{i}_template.npz');actual=sparse.load_npz(folder/f'{i}_actual.npz');libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
  parent_count,parent_diag=emit(template,np.arange(len(libs)),effects[ti],own,cell['seed']+ti+offset,'residual10',(libs,cp,bulk),atlas)
  basepr,basej=profile_jackknife_chunked(parent_count);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
  for arm in arms:
   path=dest/f'{name}__{i}__{arm}.json';arr=path.with_suffix('.npz')
   if not path.exists():
    start=time.time()
    if arm=='v66':count=parent_count.copy();diag=dict(CPM_L1=0.,fallback=False,changed_fraction=0.,bulk_exact=True,row_libraries_exact=True)
    else:
     shuffled=arm.endswith('_shuffled');kind=arm.removesuffix('_shuffled');delta=fake[ti] if shuffled else genuine[ti]
     count,diag=m.transform(parent_count,query,delta,kind,cell['seed']+ti+offset+7800000,shuffled,panel)
    np.testing.assert_array_equal(np.asarray(count.sum(0)).ravel(),np.asarray(parent_count.sum(0)).ravel());np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
    pr,pj=profile_jackknife_chunked(count);np.testing.assert_array_equal(pr,basepr)
    dd=de(count,reference,genes);metrics=direction_metrics(dd,truth,genes,own);rank_checks=spot_ranks(count,reference,dd)
    if arm=='v66':
     old=m.ROOT/'experiments/v66_assay_strength/analysis/counts'/f'{name}__{i}__{offset}__assay75.npz'
     with np.load(old,allow_pickle=False) as z:
      for field,v in [('pred',pr),('jack',pj),('pred_q',dd['q']),('pred_lfc',dd['lfc']),('truth',tp)]:np.testing.assert_allclose(z[field],v,rtol=0,atol=1e-12)
    np.savez_compressed(arr,pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
    r=dict(family=name,study=m.key(name),target=target,case_index=i,arm=arm,seed_offset=offset,cells=len(libs),exact_V66_bulk_profile=True,exact_V66_gene_totals=True,exact_original_libraries=True,manual_rank_coordinates_checked=rank_checks,baseline_full_replay_verified=arm=='v66',**metrics,**diag,arrays_sha256=sha(arr),seconds=time.time()-start)
    write(path,r);del count;gc.collect()
   r=read(path);assert sha(arr)==r['arrays_sha256']
   with np.load(arr,allow_pickle=False) as z:pr=z['pred'].copy();pj=z['jack'].copy();np.testing.assert_array_equal(pr,basepr)
   outputs[arm].append((pr,pj,tp,tj,r));print('V78',stage,offset,name,i,arm,round(r['seconds'],1),'fallback',r['fallback'],flush=True)
  del parent_count;gc.collect()
 for arm,values in outputs.items():
  pp,jp,tt,jt=(np.asarray([v[k] for v in values]) for k in range(4));rows=[v[4] for v in values]
  result=panel_metrics(pp,tt,cc,jp,jt,cj,genes,[r['target'] for r in rows],targets)
  for field in ('fidelity','reach','jaccard','nmae'):
   vals=[r[field] for r in rows if r.get(field) is not None and np.isfinite(r[field])];result[field]=float(np.mean(vals)) if vals else None
  result['fallback_fraction']=float(np.mean([r['fallback'] for r in rows]));result['CPM_L1']=float(np.mean([r['CPM_L1'] for r in rows]));result['changed_fraction']=float(np.mean([r['changed_fraction'] for r in rows]))
  write(dest/f'{name}__{arm}__panel.json',dict(family=name,study=m.key(name),arm=arm,**result))
 for arm in arms:
  panelrecord=read(dest/f'{name}__{arm}__panel.json');baseline=read(dest/f'{name}__v66__panel.json');assert panelrecord['pds']==baseline['pds']

def assess(stage,offset):
 if stage=='development':tested=m.ARMS
 else:tested=(read(m.A/'selection.json')['selected'],)
 arms=('v66',)+tuple(tested)+tuple(a+'_shuffled' for a in tested);metrics=('pds','fidelity','reach','jaccard','nmae','CPM_L1','changed_fraction','fallback_fraction','common_expression_numerator')
 studies={};summary={};dest=m.A/stage/str(offset)
 for arm in arms:
  records=[read(dest/f"{ctx['family']}__{arm}__panel.json") for ctx in m.contexts()];studies[arm]={}
  for s in ('H1','Orion','CD4'):
   studies[arm][s]={}
   for field in metrics:
    v=[r[field] for r in records if r['study']==s and r.get(field) is not None and np.isfinite(r[field])];studies[arm][s][field]=float(np.mean(v)) if v else None
  summary[arm]={}
  for field in metrics:
   v=[r[field] for r in studies[arm].values() if r[field] is not None];summary[arm][field]=float(np.mean(v)) if len(v)==3 else None
  available=[s for s in studies[arm] if studies[arm][s]['nmae'] is not None];summary[arm]['nmae_available_studies']=available;summary[arm]['nmae_available_mean']=float(np.mean([studies[arm][s]['nmae'] for s in available])) if available else None
 baseline=summary['v66'];checks={}
 for arm in tested:
  v=summary[arm];sh=summary[arm+'_shuffled'];nmae=v['nmae_available_mean'];bn=baseline['nmae_available_mean']
  checks[arm]=dict(profile_identity=abs(v['pds']-baseline['pds'])<1e-15,CPM_bound=max(r['CPM_L1'] for r in studies[arm].values())<=.0005,fidelity=v['fidelity']>=baseline['fidelity']+.002,jaccard=v['jaccard']>=baseline['jaccard']+.0005,reach=v['reach']>=baseline['reach']-.005,nmae_available=nmae is not None and bn is not None and nmae<=1.01*bn,beats_shuffle=v['fidelity']>sh['fidelity'] and v['jaccard']>=sh['jaccard'])
 eligible=[a for a in tested if all(checks[a].values())]
 selected=max(eligible,key=lambda a:(summary[a]['fidelity'],summary[a]['jaccard'])) if eligible else None
 result=dict(stage=stage,seed_offset=offset,summary=summary,studies=studies,checks=checks,eligible=eligible,selected=selected,complete_six_metric_support=False,scope='Exact raw profile/PDS identity; DE and cell-normalized means may change; partial NMAE guard uses only stated available studies. Reused public biology.')
 write(m.A/f'{stage}_{offset}_assessment.json',result)
 if stage=='development':
  assert not (m.A/'selection.json').exists();write(m.A/'selection.json',dict(selected=selected,eligible=eligible,development_sha256=sha(m.A/f'{stage}_{offset}_assessment.json')))
 print(json.dumps(dict(stage=stage,selected=selected,summary=summary,checks=checks)),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['counts','assess']);p.add_argument('--stage',choices=['development','confirmation'],required=True);p.add_argument('--source');p.add_argument('--offset',type=int,default=0);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.action=='counts':counts(args.source,args.stage,args.offset)
  else:assess(args.stage,args.offset)
