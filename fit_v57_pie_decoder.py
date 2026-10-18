"""Fit source-only three-head calibration; evaluate every eligible CD4 target untouched."""
from functools import lru_cache
import argparse,gc
import numpy as np
from sklearn.linear_model import Ridge
from v57_pie import *
from build_v57_pie_evidence import blocks

KINDS=('joint','shuffled','evidence')
@lru_cache(5)
def inputs(source):
 p=CACHE/'inputs'/(source+'.npz');assert sha(p)==read(A/'inputs'/(source+'.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
@lru_cache(5)
def predictions(source):
 p=CACHE/'predictions'/(source+'.npz');assert sha(p)==read(A/'predictions'/(source+'.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def features(source,ti,pos,kind):
 z=inputs(source);pr=predictions(source);pp=plan();j=pp['shuffle'][ti] if kind=='shuffled' else ti
 assert pr['done'][j];col=np.searchsorted(pr['positions'],pos);assert np.array_equal(pr['positions'][col],pos)
 b=z['baseline'][ti,pos];c=z['ctrl_mean'][pos];q=np.log1p(z['query_cpm'][pos]);base=[b,c,q,b*c,b*b]
 ev,_=blocks(str(z['targets'][j]),source,z['genes'][pos].tolist());ed=ev['evidence_dp'];el=ev['evidence_lfc'];ee=ev['evidence_de']
 more=[ed[:,0],ed[:,1],ed[:,2],ed[:,3],el[:,0],el[:,1],ee[:,0]-ee[:,1],ee[:,0]+ee[:,1],ee[:,2],ed[:,0]*b]
 if kind=='evidence':return np.column_stack(base+more).astype(np.float32)
 lf=np.clip(pr['lfc'][j,col],-3,3);p=pr['p_de'][j,col];dp=pr['delta_p'][j,col]
 new=[p,lf,dp,p*lf,p*dp,abs(dp),dp/(c+.1),dp*c,lf*b,dp*b,p*ed[:,0],(2*p-1)*np.sign(lf)]
 return np.column_stack(base+more+new).astype(np.float32)

def normalizer(x):
 mean=x.mean(0);sd=x.std(0);sd=np.where(sd>1e-8,np.maximum(sd,.01),np.inf)
 return mean,sd

def train(kind):
 pp=plan();xx=[];yy=[];ww=[];validation=[];groups=[]
 for source in pp['config']['calibration_sources']:
  z=inputs(source);ncase=len(pp['cohorts'][source]['targets'])
  for target in pp['cohorts'][source]['targets']:
   ti=list(z['targets']).index(target);valid=z['valid'][ti]&z['support']&np.isfinite(z['delta'][ti]).all(0)
   valid[list(z['genes']).index(target)]=False
   pos=np.flatnonzero(valid);rng=np.random.default_rng(SEED+ti);pos=rng.choice(pos,min(len(pos),768),replace=False)
   x=features(source,ti,pos,kind);y=np.column_stack([z['mean'][ti].mean(0)[pos]-z['baseline'][ti,pos],z['rank'][ti].mean(0)[pos],z['delta'][ti].mean(0)[pos]])
   good=np.isfinite(x).all(1)&np.isfinite(y).all(1);x=x[good];y=y[good]
   xx.append(x);yy.append(y);ww.extend([1/ncase]*len(x));validation.extend([split(target)==0]*len(x));groups.append({'source':source,'target':target,'observations':len(x),'validation':split(target)==0})
 x=np.concatenate(xx).astype(np.float64);y=np.concatenate(yy).astype(np.float64);w=np.array(ww);w/=w.mean();v=np.array(validation);assert v.any() and (~v).any()
 mean,sd=normalizer(x[~v]);scale=np.maximum(y[~v].std(0),[.05,.01,.01]);xt=np.clip((x-mean)/sd,-8,8);scores=[]
 for alpha in pp['config']['ridge_alphas']:
  model=Ridge(alpha=alpha).fit(xt[~v],y[~v]/scale,sample_weight=w[~v]);err=model.predict(xt[v])*scale-y[v];loss=np.average((err/scale)**2,axis=0,weights=w[v]);scores.append({'alpha':alpha,'normalized_MSE':loss.tolist(),'objective':float(loss.mean())})
 best=min(scores,key=lambda a:a['objective'])['alpha'];mean,sd=normalizer(x);scale=np.maximum(y.std(0),[.05,.01,.01]);model=Ridge(alpha=best).fit(np.clip((x-mean)/sd,-8,8),y/scale,sample_weight=w)
 folder=CACHE/'calibrators';folder.mkdir(exist_ok=True);path=folder/(kind+'.npz');np.savez_compressed(path,mean=mean,sd=sd,scale=scale,coef=model.coef_,intercept=model.intercept_)
 write(A/'calibrators'/(kind+'.json'),{'sha256':sha(path),'plan_sha256':sha(A/'plan.json'),'kind':kind,'alpha':best,'source_target_validation':scores,'groups':groups,'source_only':True,'CD4_labels_used_in_fit':False,'target_scales':scale.tolist(),'inactive_features':np.flatnonzero(~np.isfinite(sd)).tolist(),'rows':len(x)})
 print('Calibration fitted',kind,'alpha',best,'rows',len(x),flush=True)

@lru_cache(3)
def model(kind):
 path=CACHE/'calibrators'/(kind+'.npz');assert sha(path)==read(A/'calibrators'/(kind+'.json'))['sha256']
 with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def predict(source,ti,kind='joint'):
 z=inputs(source);pos=np.flatnonzero(z['support']);x=features(source,ti,pos,kind);m=model(kind)
 raw=(np.clip((x-m['mean'])/m['sd'],-8,8)@m['coef'].T+m['intercept'])*m['scale']
 mean=z['baseline'][ti].copy();mean[pos]+=np.clip(raw[:,0],-1,1);mean=np.clip(mean,-3,3)
 rank=np.zeros_like(mean);rank[pos]=np.clip(raw[:,1],-.25,.25)
 delta=np.full_like(mean,np.nan);delta[pos]=raw[:,2]
 own=list(z['genes']).index(str(z['targets'][ti]));mean[own]=0;rank[own]=0
 return mean,rank,delta

def evaluate():
 import pandas as pd
 pp=plan();rows=[]
 for source in SOURCES:
  z=inputs(source);pr=predictions(source)
  for target in pp['cohorts'][source]['targets']:
   ti=list(z['targets']).index(target);mask=z['valid'][ti]&z['support']&np.isfinite(z['delta'][ti]).all(0);pos=np.flatnonzero(mask)
   ym=z['mean'][ti].mean(0);yr=z['rank'][ti].mean(0);yd=z['delta'][ti].mean(0)
   cols=np.searchsorted(pr['positions'],pos)
   row={'source':source,'target':target,'seen':target in pp['seen_targets'],'transfer':source.startswith('CD4'),'genes':len(pos),'baseline_mean_MSE':float(np.mean((z['baseline'][ti,pos]-ym[pos])**2)),'zero_rank_MSE':float(np.mean(yr[pos]**2)),'zero_delta_MSE':float(np.mean(yd[pos]**2)),'raw_PIE_delta_MSE':float(np.mean((pr['delta_p'][ti,cols]-yd[pos])**2)),'mean_view_disagreement':float(np.mean((z['mean'][ti,0,pos]-z['mean'][ti,1,pos])**2)),'rank_view_disagreement':float(np.mean((z['rank'][ti,0,pos]-z['rank'][ti,1,pos])**2))}
   for kind in KINDS:
    mu,rank,delta=predict(source,ti,kind);row.update({kind+'_mean_MSE':float(np.mean((mu[pos]-ym[pos])**2)),kind+'_rank_MSE':float(np.mean((rank[pos]-yr[pos])**2)),kind+'_delta_MSE':float(np.mean((delta[pos]-yd[pos])**2))})
   rows.append(row)
  print('Response evaluated',source,flush=True)
 frame=pd.DataFrame(rows);frame.to_csv(A/'response_cases.csv',index=False)
 frame.groupby(['source','seen']).mean(numeric_only=True).to_csv(A/'response_summary.csv')
 write(A/'response_evaluation.json',{'complete':True,'cases':len(rows),'held_CD4_cases':int(frame.transfer.sum()),'response_cases_sha256':sha(A/'response_cases.csv'),'summary_sha256':sha(A/'response_summary.csv'),'plan_sha256':sha(A/'plan.json'),'interpretation':'Equal target metrics, disjoint-view labels. Orion/H1 overlap pretrained model training and are calibration diagnostics. CD4 excluded from decoder fit. No count or private leaderboard gain implied.'})
 print(frame.groupby(['source','seen'])[['baseline_mean_MSE','joint_mean_MSE','zero_rank_MSE','joint_rank_MSE','zero_delta_MSE','raw_PIE_delta_MSE','joint_delta_MSE']].mean().to_string(),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['train','evaluate']);p.add_argument('--kind',choices=KINDS);args=p.parse_args()
 if args.stage=='train':
  for k in ([args.kind] if args.kind else KINDS):train(k)
 else:evaluate()
