"""End-to-end count tests: V24, calibrated means, joint decoder and shuffled control."""
import argparse,gc,time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import ndtri
from sklearn.metrics import average_precision_score
from threadpoolctl import threadpool_limits
from v57_pie import *
from fit_v57_pie_decoder import inputs,predictions,predict
from v25_methods import raw
from v16_methods import emit as old_emit
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile,profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model


def emit_joint(template,effect,rank,own,seed,atlas):
 depth,p,b=moments(template)[:3];libs=np.asarray(depth,dtype=np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=b*factor;wb/=wb.sum()
 original=template.toarray().astype(np.float64);empirical=original*factor;empirical*=np.divide(libs,empirical.sum(1))[:,None]
 normalized=empirical/libs[:,None]*1e6;control=original/libs[:,None]*1e6
 delta=(wp-p)*1e6;theta=np.clip(abs(rank),.005,.249)
 desired=np.maximum((delta/ndtri(.5+theta))**2-control.var(0),0)
 variance=normalized.var(0);poisson_variance=(.65**2)*1e12*wp*np.mean(1/libs)
 weight=np.sqrt(np.maximum(desired-poisson_variance,0)/np.maximum(variance,1e-12));weight=np.clip(weight,0,1)
 active=(abs(rank)>.005)&(abs(effect)>1e-6)&(np.sign(rank)==np.sign(effect))
 weight=np.where(active,weight,1.)
 if own is not None:weight[own]=1.
 expected=libs[:,None]*wp;pois=np.random.default_rng(seed).poisson(expected)
 beta=.65*(1-weight);initial=(1-weight-beta)[None,:]*expected+beta[None,:]*pois+weight[None,:]*empirical
 assert (initial>=-1e-8).all();initial=np.maximum(initial,0);initial/=initial.sum(1,keepdims=True)
 count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
 assert (count>=0).all() and np.issubdtype(count.dtype,np.integer);np.testing.assert_array_equal(count.sum(1),libs)
 return sparse.csr_matrix(count),{'active_rank_coordinates':int(active.sum()),'empirical_weight_mean':float(weight.mean()),'empirical_weight_active_mean':float(weight[active].mean()) if active.any() else None,'mean_probability_L1':float(abs((count/libs[:,None]).mean(0)-wp).sum()),'bulk_probability_L1':float(abs(count.sum(0)/libs.sum()-wb).sum())}

def digest(x):
 import hashlib
 x=sparse.csr_matrix(x);h=hashlib.sha256()
 for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
 return h.hexdigest()

def run(source):
 pp=plan();z=inputs(source);genes=z['genes'];x,meta,measured,receipt=raw(source,genes);labels=meta.gene_target.to_numpy(str);ntc=np.flatnonzero(labels=='Non-Targeting');atlas,_=load_atlas_model()
 folder=A/'counts';folder.mkdir(exist_ok=True);artifact=CACHE/'count_statistics';artifact.mkdir(exist_ok=True)
 query_pred=predictions(source);alltargets=list(z['targets'])
 for target in pp['cohorts'][source]['count_pilot']:
  ti=alltargets.index(target);own=list(genes).index(target);base=z['baseline'][ti];mu,rank,_=predict(source,ti,'joint');smu,srank,_=predict(source,ti,'shuffled');arms={'v24':(base,None),'pie_mean':(mu,None),'pie_joint':(mu,rank),'shuffled_joint':(smu,srank)}
  for seed in (pp['config']['count_seeds'] if source.startswith('CD4') else pp['config']['count_seeds'][:1]):
   path=folder/f'{source}_{target}_{seed}.json'
   if path.exists():assert read(path)['plan_sha256']==sha(A/'plan.json');continue
   rng=np.random.default_rng(seed+ti);order=rng.permutation(ntc);template_ids=order[:400];reference_ids=order[400:2400];tr=np.flatnonzero(labels==target);actual_ids=rng.permutation(tr)[:400]
   assert len(actual_ids)>=16 and len(template_ids)==400 and not set(template_ids)&set(reference_ids)
   template=x[template_ids];reference=x[reference_ids];actual=x[actual_ids];truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual);cp,cj=profile_jackknife_chunked(reference)
   prepared=moments(template)[:3];records=[];arrays={'truth_profile':tp,'truth_jackknife':tj,'control_profile':cp,'control_jackknife':cj,'truth_q':truth['q'],'truth_lfc':truth['lfc']}
   pos=np.flatnonzero(z['support']&(np.arange(len(genes))!=own));col=np.searchsorted(query_pred['positions'],pos);truth_de=(truth['q'][pos]<.05)&truth['eligible'][pos];prob=query_pred['p_de'][ti,col]
   headap=float(average_precision_score(truth_de,prob)) if truth_de.any() else None
   for arm,(effect,rk) in arms.items():
    begin=time.time()
    if rk is None:count,diag=old_emit(template,np.arange(400),effect,own,seed+ti,'residual10',prepared=prepared,atlas=atlas)
    else:count,diag=emit_joint(template,effect,rk,own,seed+ti,atlas)
    stat=de(count,reference,genes);score=direction_metrics(stat,truth,genes,own);pred,pj=profile_jackknife_chunked(count)
    axis=~np.isin(genes,z['targets']);mse=float(np.mean((pred[axis]-tp[axis])**2));unbiased=float(np.mean(((pred-tp)**2-pj-tj)[axis]))
    records.append({'source':source,'target':target,'seen':target in pp['seen_targets'],'transfer':source.startswith('CD4'),'seed':seed,'arm':arm,'actual_cells':len(actual_ids),'generated_cells':400,'reference_cells':len(reference_ids),'profile_MSE':mse,'profile_error_unbiased':unbiased,'count_sha256':digest(count),'raw_PIE_pDE_average_precision':headap,**score,**diag,'seconds':time.time()-begin})
    arrays[arm+'_profile']=pred;arrays[arm+'_jackknife']=pj;arrays[arm+'_q']=stat['q'];arrays[arm+'_lfc']=stat['lfc']
    print('Counts',source,target,seed,arm,'FID',round(score['fidelity'],4),'MSE',round(mse,5),round(time.time()-begin,1),'sec',flush=True)
    del count;gc.collect()
   file=artifact/path.with_suffix('.npz').name;np.savez_compressed(file,**arrays)
   write(path,{'plan_sha256':sha(A/'plan.json'),'records':records,'statistics_sha256':sha(file),'template_rows':template_ids.tolist(),'reference_rows':reference_ids.tolist(),'actual_rows':actual_ids.tolist(),'source_record':receipt,'interpretation':'Source overlapping studies are calibration diagnostics. CD4 not used in fitting. Exact integer counts scored after row normalization/rounding; no official anchors or private-score claim.'})
 # Real NTC identity anchor tested without a perturbation or a predicted rank effect.
 nullpath=folder/(source+'_null.json')
 if not nullpath.exists():
  rng=np.random.default_rng(SEED);order=rng.permutation(ntc);template=x[order[:400]];ref=x[order[400:2400]];effect=np.zeros(len(genes));rr=np.zeros(len(genes));nn=[]
  for arm in ['v24','pie_joint']:
   if arm=='v24':count,diag=old_emit(template,np.arange(400),effect,None,SEED,'residual10',prepared=moments(template)[:3],atlas=atlas)
   else:count,diag=emit_joint(template,effect,rr,None,SEED,atlas)
   stat=de(count,ref,genes);nn.append({'arm':arm,'DE_calls':int((stat['q']<.05).sum()),'count_equals_real_template':bool((count!=template).nnz==0),**diag})
  write(nullpath,{'plan_sha256':sha(A/'plan.json'),'source':source,'records':nn,'reference_cells':len(ref),'identity_anchor_scope':'Zero mean change and zero rank effect only; no inference about intervention accuracy'})


def summarize():
 pp=plan();files=sorted((A/'counts').glob('*.json'));records=[]
 for f in files:
  if f.stem.endswith('_null'):continue
  r=read(f);assert r['plan_sha256']==sha(A/'plan.json');records.extend(r['records'])
 frame=pd.DataFrame(records);frame.to_csv(A/'count_cases.csv',index=False)
 frame.groupby(['source','seen','arm']).mean(numeric_only=True).to_csv(A/'count_summary.csv')
 # Paired raw-score changes; no conversion to official scaled scores without anchors.
 paired=[]
 for (source,target,seed),g in frame.groupby(['source','target','seed']):
  b=g.set_index('arm');base=b.loc['v24']
  for arm in ['pie_mean','pie_joint','shuffled_joint']:
   r=b.loc[arm];paired.append({'source':source,'target':target,'seed':int(seed),'seen':bool(r.seen),'arm':arm,**{k:float(r[k]-base[k]) if pd.notna(r[k]) and pd.notna(base[k]) else np.nan for k in ['fidelity','reach','jaccard','nmae','profile_MSE','profile_error_unbiased']}})
 pd.DataFrame(paired).to_csv(A/'count_paired_changes.csv',index=False)
 expected=sum(len(v['count_pilot'])*(2 if s.startswith('CD4') else 1)*4 for s,v in pp['cohorts'].items());complete=len(frame)==expected
 write(A/'count_evaluation.json',{'complete':complete,'comparisons':len(frame),'expected':expected,'plan_sha256':sha(A/'plan.json'),'counts_sha256':sha(A/'count_cases.csv'),'summary_sha256':sha(A/'count_summary.csv'),'D3_opened':False,'official_private_score':None})
 print(frame.groupby(['source','arm'])[['fidelity','reach','jaccard','nmae','profile_MSE']].mean().to_string(),flush=True)

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['run','summarize']);q.add_argument('--source',choices=SOURCES);args=q.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='summarize':summarize()
  else:
   for s in ([args.source] if args.source else SOURCES):run(s)
