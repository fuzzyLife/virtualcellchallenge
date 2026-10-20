"""V70 paired confirmation counts with V66 and matched shuffled evidence."""
import argparse,gc,json,time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v76_response_dilution as m
from v57_pie import read,write,sha
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics


def counts(name):
 m.check();selected=read(m.A/'selection.json')['selected'];assert selected
 ctx=next(c for c in m.contexts() if c['family']==name);folder,_=m.checked_folder(name)
 d=m.data();genes=d['genes'];targets=d['targets'];dest=m.A/'counts_repeat';dest.mkdir(exist_ok=True)
 path=m.A/'predictions'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256']
 with np.load(path,allow_pickle=False) as z:pred={a:z[a].copy() for a in ('v66',selected,selected+'_shuffled')}
 with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy();bulk=z['bulk'].copy()
 reference=sparse.load_npz(folder/'reference.npz');cc,cj=profile_jackknife_chunked(reference)
 atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();arms=list(pred);outputs={a:[] for a in arms}
 for i in m.cases(ctx):
  cell=ctx['cells'][i];ti=list(targets).index(cell['target']);own=list(genes).index(cell['target'])
  template=sparse.load_npz(folder/f'{i}_template.npz');actual=sparse.load_npz(folder/f'{i}_actual.npz');libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
  truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual);axis=np.arange(len(genes))!=own
  for arm in arms:
   path=dest/f'{name}__{i}__{arm}.json';arr=path.with_suffix('.npz')
   if not path.exists():
    start=time.time()
    count,diag=emit(template,np.arange(len(libs)),pred[arm][ti],own if arm=='v66' else None,cell['seed']+ti+620000,'residual10',(libs,cp,bulk),atlas)
    np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);assert count.dtype.kind in 'iu' and count.data.min(initial=0)>=0
    dd=de(count,reference,genes);pr,pj=profile_jackknife_chunked(count);metrics=direction_metrics(dd,truth,genes,own)
    np.savez_compressed(arr,pred_p=dd['p'],truth_p=truth['p'],pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
    r=dict(family=name,study=m.study(name),target=cell['target'],case_index=i,arm=arm,cells=len(libs),**metrics,**diag)
    from audit_v76_counts import spot_ranks
    checked=spot_ranks(count,reference,dd)
    r['manual_rank_coordinates_checked']=checked
    del count;gc.collect()
    if arm=='v66':
     old=m.ROOT/'experiments/v66_assay_strength/analysis/counts'/f'{name}__{i}__620000__assay75.json';oldrec=read(old);assert sha(old.with_suffix('.npz'))==oldrec['arrays_sha256']
     with np.load(old.with_suffix('.npz'),allow_pickle=False) as z:
      for field,value in [('pred',pr),('jack',pj),('truth',tp),('truth_jack',tj),('pred_q',dd['q']),('pred_lfc',dd['lfc']),('pred_eligible',dd['eligible'])]:
       np.testing.assert_allclose(z[field],value,rtol=0,atol=1e-12)
     r['baseline_full_replay_verified']=True;r['original_baseline_receipt_sha256']=sha(old)
    r['arrays_sha256']=sha(arr);r['seconds']=time.time()-start;write(path,r)
   r=read(path);assert sha(arr)==r['arrays_sha256']
   with np.load(arr,allow_pickle=False) as z:pr=z['pred'].copy();pj=z['jack'].copy();np.testing.assert_array_equal(z['truth'],tp)
   r['heldout_mean_case']=i in m.cases(ctx)[3:]
   r['profile_square_error']=float(np.sum((pr[axis]-tp[axis])**2));outputs[arm].append((pr,pj,tp,tj,r))
   print('V70 counts',name,i,arm,flush=True)
 for arm,vals in outputs.items():
  pp,jp,tt,jt=(np.asarray([v[k] for v in vals]) for k in range(4));rows=[v[4] for v in vals]
  panel=panel_metrics(pp,tt,cc,jp,jt,cj,genes,[r['target'] for r in rows],targets)
  for metric in ('fidelity','reach','jaccard','nmae','profile_square_error'):
   vv=[r[metric] for r in rows if r['heldout_mean_case'] and r.get(metric) is not None and np.isfinite(r[metric])];panel[metric]=float(np.mean(vv)) if vv else None
  write(dest/f'{name}__{arm}__panel.json',dict(family=name,study=m.study(name),arm=arm,PDS_unique_targets=len(vals),profile_DE_targets=len(m.cases(ctx)[3:]),**panel))


def assess():
 selected=read(m.A/'selection.json')['selected'];arms=('v66',selected,selected+'_shuffled');metrics=('pds','fidelity','reach','jaccard','nmae','profile_square_error');studies={}
 for arm in arms:
  records=[read(m.A/'counts_repeat'/f"{c['family']}__{arm}__panel.json") for c in m.contexts()]
  studies[arm]={}
  for s in ('H1','Orion','CD4'):
   studies[arm][s]={}
   for metric in metrics:
    v=[r[metric] for r in records if r['study']==s and r.get(metric) is not None and np.isfinite(r[metric])]
    studies[arm][s][metric]=float(np.mean(v)) if v else None
 summary={a:{k:float(np.mean([s[k] for s in studies[a].values()])) if all(s[k] is not None for s in studies[a].values()) else None for k in metrics} for a in arms}
 ratios={s:studies[selected][s]['profile_square_error']/studies['v66'][s]['profile_square_error'] for s in studies[selected]}
 base=summary['v66'];p=summary[selected];shuffle=summary[selected+'_shuffled'];meanratio=float(np.mean(list(ratios.values())))
 checks=dict(profile=meanratio<=.995,worst_study=max(ratios.values())<=1.02,
  pds=p['pds']>=base['pds']-.005,fidelity=p['fidelity']>=base['fidelity']-.002,
  reach=p['reach']>=base['reach']-.005,jaccard=p['jaccard']>=base['jaccard']-.0005,
  nmae=p['nmae'] is not None and base['nmae'] is not None and p['nmae']<=1.01*base['nmae'],
  beats_shuffle=p['profile_square_error']<shuffle['profile_square_error'])
 result=dict(seed_offset=620000,selected=selected,checks=checks,release_supported=all(checks.values()),summary=summary,studies=studies,profile_ratios=ratios,paired_profile_ratio=meanratio,missing_metrics=[k for k in metrics if p[k] is None],D3_opened=False,scope='Reused public biological cohorts; fixed remaining unique targets, one technical seed; no private score forecast.')
 write(m.A/'count_assessment_repeat.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['counts','assess']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='counts':counts(a.source)
  else:assess()
