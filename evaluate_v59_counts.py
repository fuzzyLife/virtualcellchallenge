"""V59 paired count comparisons; no upload, no private-score estimate."""
import argparse,gc,time
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from v59_state_composition import *
from v25_methods import raw
from v16_methods import emit
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
COUNT_ARMS=(*ARMS,'state_cells50')

def state_counts(x,pool,state,tilt,within,base,own,seed):
 rng=np.random.default_rng(seed);pool=np.asarray(pool)
 if own is None and not np.any(tilt) and not np.any(within) and not np.any(base):
  rows=rng.choice(pool,400,replace=False);return x[rows].copy(),rows
 first=rng.choice(pool,200,replace=False);remain=np.setdiff1d(pool,first);w=np.exp(tilt[state[remain]]);key=rng.exponential(1/w);second=remain[np.argsort(key,kind='stable')[:200]];rows=np.concatenate([first,second]);assert len(set(rows))==400
 template=x[rows].toarray().astype(np.int64);libs=template.sum(1)
 effect=np.vstack([np.tile(base,(200,1)),within[state[second]]]);factor=np.exp2(np.clip(effect,-3,3));
 if own is not None:factor[:,own]=.1
 expected=template*factor;expected*=np.divide(libs,expected.sum(1))[:,None]
 rate=np.divide(expected,template,out=np.ones_like(expected),where=template>0);removed=rng.binomial(template,1-np.minimum(rate,1));count=template-removed;extra=np.maximum(expected-template,0)
 for i,n in enumerate(removed.sum(1)):
  if n:
   assert extra[i].sum()>0;count[i]+=rng.multinomial(int(n),extra[i]/extra[i].sum())
 assert (count>=0).all();np.testing.assert_array_equal(count.sum(1),libs)
 return sparse.csr_matrix(count),rows

def run(source):
 pp=plan();genes,targets=axes();predpath=CACHE/'predictions'/(source+'.npz');assert sha(predpath)==read(A/'predictions'/(source+'.json'))['sha256'];pred=load(predpath)
 x,meta,measured,sourcerec=raw(source,genes);label=meta.gene_target.to_numpy(str);ntc=np.flatnonzero(label=='Non-Targeting');atlas,_=load_atlas_model();mech=mechanism(GROUP[source],source);state=np.zeros(x.shape[0],int);state[mech['rows']]=mech['state']
 records_dir=A/'counts';records_dir.mkdir(exist_ok=True);statsdir=CACHE/'count_statistics';statsdir.mkdir(exist_ok=True)
 for seed in (CONFIG['count_seeds'] if source.startswith('CD4') else CONFIG['count_seeds'][:1]):
  order=np.random.default_rng(seed).permutation(ntc);poolsize=min(max(400,len(order)//2),len(order)-400);pool=order[:poolsize];reference_ids=order[poolsize:poolsize+2000];template_ids=np.random.default_rng(seed+1).choice(pool,400,replace=False)
  assert len(reference_ids)>=400 and not set(pool)&set(reference_ids)
  reference=x[reference_ids];template=x[template_ids];prepared=moments(template)[:3];cp,cj=profile_jackknife_chunked(reference)
  for target in pp['cohorts'][source]['count_pilot']:
   ti=list(targets).index(target);own=list(genes).index(target);path=records_dir/f'{source}_{target}_{seed}.json'
   if path.exists():assert read(path)['plan_sha256']==sha(A/'plan.json');continue
   actual_ids=np.random.default_rng(seed+ti).permutation(np.flatnonzero(label==target))[:400];actual=x[actual_ids];truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual);records=[];arrays={'truth_profile':tp,'truth_jackknife':tj,'control_profile':cp,'control_jackknife':cj,'truth_q':truth['q'],'truth_lfc':truth['lfc']}
   for arm in COUNT_ARMS:
    begin=time.time()
    if arm=='state_cells50':
     count,rows=state_counts(x,pool,state,pred['tilt'][ti],pred['within'][ti],pred['v24'][ti],own,seed+ti);diag={'actual_control_rows':rows.tolist()}
    else:count,diag=emit(template,np.arange(400),pred[arm][ti],own,seed+ti,'residual10',prepared,atlas)
    stat=de(count,reference,genes);score=direction_metrics(stat,truth,genes,own);profile,pj=profile_jackknife_chunked(count);axis=np.arange(len(genes))!=own;nonpanel=~np.isin(genes,targets)
    record=dict(source=source,study=GROUP[source],target=target,seed=seed,arm=arm,actual_cells=len(actual_ids),reference_cells=len(reference_ids),profile_MSE=float(np.mean((profile[axis]-tp[axis])**2)),nonpanel_profile_MSE=float(np.mean((profile[nonpanel]-tp[nonpanel])**2)),profile_error_unbiased=float(np.mean(((profile-tp)**2-pj-tj)[axis])),**score,seconds=time.time()-begin)
    records.append(record);arrays[arm+'_profile']=profile;arrays[arm+'_jackknife']=pj;arrays[arm+'_q']=stat['q'];arrays[arm+'_lfc']=stat['lfc'];print('Counts',source,target,seed,arm,'FID',round(score['fidelity'],4),'MSE',round(record['profile_MSE'],5),round(time.time()-begin,1),'sec',flush=True);del count;gc.collect()
   artifact=statsdir/path.with_suffix('.npz').name;save(artifact,**arrays)
   write(path,dict(plan_sha256=sha(A/'plan.json'),records=records,statistics_sha256=sha(artifact),template_rows=template_ids.tolist(),template_pool=pool.tolist(),reference_rows=reference_ids.tolist(),actual_rows=actual_ids.tolist(),state_cell_template_rows=diag.get('actual_control_rows'),source_record=sourcerec,D3_opened=False))
 null=A/'counts'/(source+'_null.json')
 if not null.exists():
  rng=np.random.default_rng(SEED);order=rng.permutation(ntc);t=x[order[:400]];r=x[order[400:2400]];zero=np.zeros(len(genes));count,_=emit(t,np.arange(400),zero,None,SEED,'residual10',moments(t)[:3],atlas);s=de(count,r,genes);real=de(t,r,genes)
  write(null,dict(source=source,residual10_DE=int((s['q']<.05).sum()),real_control_DE=int((real['q']<.05).sum()),state_zero_anchor='At zero intervention return unedited real control rows; no artificial knockdown applied',plan_sha256=sha(A/'plan.json')))

def summarize():
 pp=plan();records=[]
 for path in sorted((A/'counts').glob('*.json')):
  if path.stem.endswith('_null'):continue
  r=read(path);assert r['plan_sha256']==sha(A/'plan.json');records.extend(r['records'])
 df=pd.DataFrame(records);df.to_csv(A/'count_cases.csv',index=False);summary=df.groupby(['study','source','arm']).mean(numeric_only=True);summary.to_csv(A/'count_summary.csv');metrics=['fidelity','reach','jaccard','nmae','profile_MSE'];print(summary[metrics].to_string(),flush=True)
 expected=sum(len(pp['cohorts'][s]['count_pilot'])*(2 if s.startswith('CD4') else 1)*len(COUNT_ARMS) for s in SOURCES);assert len(df)==expected
 wide=df.groupby(['study','source','arm'])[metrics].mean().groupby(['study','arm']).mean();checks={}
 for study in OUTERS:
  z=wide.loc[study];b=z.loc['v24'];j=z.loc['joint50'];sh=z.loc['shuffled_joint50'];checks[study]=dict(expression_better=bool(j.profile_MSE<b.profile_MSE),DE_gains=int(sum(j[k]>b[k] for k in ['fidelity','reach','jaccard'])),FID_guard=bool(j.fidelity>=b.fidelity-.005),JAC_guard=bool(j.jaccard>=b.jaccard-.005),beats_shuffle_FID=bool(j.fidelity>sh.fidelity),beats_shuffle_expression=bool(j.profile_MSE<sh.profile_MSE))
 count_pass=all(v['expression_better'] and v['DE_gains']>=2 and v['FID_guard'] and v['JAC_guard'] and v['beats_shuffle_FID'] and v['beats_shuffle_expression'] for v in checks.values());response=read(A/'response_assessment.json')
 write(A/'count_assessment.json',dict(complete=True,comparisons=len(df),expected=expected,checks=checks,count_pass=bool(count_pass),response_pass=response['passed'],release_supported=bool(count_pass and response['passed']),count_sha256=sha(A/'count_cases.csv'),plan_sha256=sha(A/'plan.json'),D3_opened=False,private_score=None))
if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['run','summarize']);q.add_argument('--source',choices=SOURCES);args=q.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='run':run(args.source)
  else:summarize()

