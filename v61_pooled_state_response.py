"""V61: pool matched state sufficient statistics BEFORE nonlinear shrinkage.

V59/V60 methods/results remain untouched. The state basis is identical to V60;
this version changes the statistical unit for abundance/conditional estimation.
"""
import argparse,gc,shutil
from functools import lru_cache
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import v60_depth_adjusted_states as parent
m=parent.m
A=m.ROOT/'experiments/v61_pooled_state_response/analysis'
CACHE=m.Path(m.tempfile.gettempdir())/'vcc_v61_pooled_state_response'


def aggregate(lanes):
 """Lane tuples (nt,nc,tm,cm,tm0,cm0); preserve target lane matching."""
 nt=np.sum([v[0] for v in lanes],axis=0);total=nt.sum();ng=lanes[0][2].shape[1]
 ts=np.zeros((4,ng));cs=np.zeros_like(ts);expected=np.zeros(4);control_sums=np.zeros_like(ts);actual=np.zeros(ng);ctrl=np.zeros(ng);ncontrol=0.
 for n,c,t,cm,tm0,cm0 in lanes:
  weight=n.sum();expected_lane=weight*c/c.sum();ts+=n[:,None]*t;cs+=n[:,None]*cm;expected+=expected_lane;control_sums+=expected_lane[:,None]*cm;actual+=weight*tm0;ctrl+=weight*cm0;ncontrol+=c.sum()
 actual/=total;ctrl/=total
 cpool=np.divide(control_sums,expected[:,None],out=np.tile(ctrl,(4,1)),where=expected[:,None]>0)
 treatment=np.divide(ts,nt[:,None],out=np.tile(actual,(4,1)),where=nt[:,None]>0)
 matched=np.divide(cs,nt[:,None],out=cpool.copy(),where=nt[:,None]>0)
 nc=expected/expected.sum()*ncontrol
 tilt,within,valid=m.estimate(nt,nc,treatment,matched)
 return dict(nt=nt,nc=nc,cm=cpool,actual=actual,control=ctrl,tilt=tilt,within=within,valid=valid,n=total)


def mechanisms(outer,source):
 m.plan();path=CACHE/'mechanisms'/outer/(source+'.npz');receipt=A/'mechanisms'/outer/(source+'.json')
 if receipt.exists():assert m.sha(path)==m.read(receipt)['sha256'];return
 genes,targets=m.axes();x,meta=m.source_data(source);state,trace=m.assign(x,meta,m.state_model(outer));lookup=np.full(int(meta['rows'].max())+1,-1,int);lookup[meta['rows']]=np.arange(len(meta['rows']))
 within=np.zeros((300,4,len(genes)),np.float32);valid=np.zeros_like(within,bool);tilt=np.zeros((300,4),np.float32);composition=np.zeros((300,len(genes)),np.float32);ncase=np.zeros(300,int);oracle=[];nonpanel=~np.isin(genes,targets)
 for group in m.source_meta(source)['groups']:
  ti=group['target_index'];target=group['target'];own=list(genes).index(target);views=[]
  for view in range(2):
   lanes=[]
   for lane in group['lanes']:
    tr=lookup[np.asarray(lane['treated_views'][view],int)];cr=lookup[np.asarray(lane['control_views'][view],int)]
    assert (tr>=0).all() and (cr>=0).all() and not set(tr)&set(cr);assert np.all(meta['label'][tr]==target) and np.all(meta['label'][cr]=='Non-Targeting')
    nt,tm,t0=m.group_summary(x,state,tr);nc,cm,c0=m.group_summary(x,state,cr);lanes.append((nt,nc,tm,cm,t0,c0))
   v=aggregate(lanes);v['valid'][:,~meta['measured']]=False;v['valid'][:,own]=False;v['within'][~v['valid']]=0;views.append(v)
  den=sum(v['valid']*v['nt'][:,None] for v in views);num=sum(v['within']*v['valid']*v['nt'][:,None] for v in views);within[ti]=np.divide(num,den,out=np.zeros_like(num),where=den>0);valid[ti]=den>0;total=sum(v['n'] for v in views);tilt[ti]=sum(v['tilt']*v['n'] for v in views)/total;composition[ti]=sum(m.logratio(m.reconstruct(v['nc'],v['cm'],v['tilt']),v['control'])*v['n'] for v in views)/total;ncase[ti]=int(total)
  for fit_view in range(2):
   train=views[fit_view];test=views[1-fit_view];axis=meta['measured']&(test['control']>5)&nonpanel;axis[own]=False;truth=m.logratio(test['actual'],test['control'])
   pred={'zero':np.zeros(len(genes)),'direct_mean':m.logratio(train['actual'],train['control']),'composition':m.logratio(m.reconstruct(test['nc'],test['cm'],train['tilt']),test['control']),'within':m.logratio(m.reconstruct(test['nc'],test['cm'],np.zeros(4),train['within']),test['control']),'joint':m.logratio(m.reconstruct(test['nc'],test['cm'],train['tilt'],train['within']),test['control']),'state_permuted':m.logratio(m.reconstruct(test['nc'],test['cm'],np.roll(train['tilt'],1)),test['control'])}
   for arm,y in pred.items():oracle.append(dict(source=source,target=target,lane='pooled matched lanes',fit_view=fit_view,arm=arm,genes=int(axis.sum()),MSE=float(np.mean((y[axis]-truth[axis])**2)),treated_cells=int(test['n']),occupancy_view_difference=float(np.linalg.norm(train['tilt']-test['tilt']))))
 ntc=np.flatnonzero(meta['label']=='Non-Targeting');nc,cm,ctrl=m.group_summary(x,state,ntc)
 digest=m.save(path,within=within,valid=valid,tilt=tilt,source_composition=composition,n=ncase,profile=ctrl,control_state_counts=nc,control_state_means=cm,state=state,rows=meta['rows'],measured=meta['measured'])
 file=A/'oracle'/outer/(source+'.csv');file.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(oracle).to_csv(file,index=False)
 m.write(receipt,dict(sha256=digest,plan_sha256=m.sha(A/'plan.json'),outer=outer,source=source,targets=int((ncase>0).sum()),oracle_sha256=m.sha(file),state_transform_control_rows=trace,training_source=outer=='all' or m.GROUP[source]!=outer,pooling='Match controls within lane/state, sum sufficient statistics across lanes, then apply density pseudocount/conditional shrinkage. Recipient outputs diagnostic only.',D3_opened=False));print('V61 pooled mechanisms',outer,source,int((ncase>0).sum()),flush=True);del x;gc.collect()


def freeze():
 assert not (A/'plan.json').exists();p=m.read(parent.A/'plan.json');p['version']=61;p['parent_plan_sha256']=m.sha(parent.A/'plan.json');p['dependencies']['v61_pooled_state_response.py']=m.sha(__file__);p['dependencies']['run_v61_experiment.py']=m.sha(m.ROOT/'run_v61_experiment.py');p['dependencies']['test_v61_pooling.py']=m.sha(m.ROOT/'test_v61_pooling.py')
 p['pooling']='Within each target/view, pool treated state counts/sums and matched-lane control state means before estimating nonlinear ratios/shrinkage. Conditional controls weight each lane by treated counts in that state. Abundance controls weight each lane by its total treated count and NTC state frequency, rescaled to actual sampled NTC total for density pseudocounts. Two view mechanisms combine with treated-state count weights. No recipient labels enter transfer.'
 p['oracle']='Fit one complete matched population view and reconstruct the other using its controls; compute response error only AFTER population aggregation. V59/V60 lane-error diagnostic was underpowered and cannot be compared numerically with this population-level error. Earlier scores and files remain intact.'
 p['reason']='Some V52 target/lane views have only two treated cells. Estimating/shrinking each lane-state effect before pooling nearly erases the response. The alternative pools the sufficient statistics while retaining batch/state matched controls.'
 p['release_candidate']='Fixed V61 joint50 primary; diagnostic arms cannot be selected post hoc. V59/V60 retained as independent comparators. A user-requested experimental upload file, if generated after failed gates, must explicitly report failure and cannot be labeled improved or overwrite V24.'
 m.write(A/'plan.json',p);print('V61 frozen',m.sha(A/'plan.json'),flush=True)


def activate():
 parent.activate();m.A=A;m.CACHE=CACHE
 @lru_cache(1)
 def plan():
  p=m.read(A/'plan.json')
  for n,h in p['dependencies'].items():assert m.sha(m.ROOT/n)==h,'V61 dependency changed: '+n
  return p
 m.plan=plan;m.mechanisms=mechanisms
 return m


def states(outer):
 source=parent.CACHE/'states'/(outer+'.npz');rfile=parent.A/'states'/(outer+'.json')
 if source.exists() and rfile.exists():
  r=m.read(rfile);assert m.sha(source)==r['sha256'];dest=CACHE/'states'/source.name;dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():assert m.sha(dest)==r['sha256']
  else:shutil.copyfile(source,dest)
  r.update(parent_state_record_sha256=m.sha(rfile),parent_state_plan_sha256=r['plan_sha256'],plan_sha256=m.sha(A/'plan.json'),basis_reused_bitwise=True);m.write(A/'states'/rfile.name,r);print('V61 reuses exact V60 state basis',outer,flush=True)
 else:m.fit_states(outer)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','states','mechanisms','transfer','evaluate','count','count_summary','plot']);p.add_argument('--source',choices=m.SOURCES);p.add_argument('--outer',choices=(*m.OUTERS,'all'));a=p.parse_args()
 if a.stage=='freeze':freeze()
 else:
  activate()
  with threadpool_limits(limits=2):
   if a.stage=='states':states(a.outer)
   elif a.stage=='mechanisms':mechanisms(a.outer,a.source)
   elif a.stage=='transfer':m.transfer(a.source)
   elif a.stage=='evaluate':m.evaluate()
   elif a.stage in ['count','count_summary']:
    import evaluate_v59_counts as counts
    if a.stage=='count':counts.run(a.source)
    else:counts.summarize()
   else:
    import plot_v59_states as plot
    plot.activate=lambda:None;plot.main()
