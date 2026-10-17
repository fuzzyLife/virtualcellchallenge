"""Train and evaluate V52 with whole-study exclusions and explicit controls."""
import argparse,copy,gc,json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import joblib,numpy as np,torch
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
import v52_joint_response as m

def loss_values(pred,y,prior,w,arm):
 den=torch.clamp(torch.mean((prior-y)**2*w,dim=0),min=torch.tensor([.04,.0004],device=y.device))
 each=torch.mean((pred-y)**2*w,dim=0)/den
 return each[0] if arm=='mean_only' else each[1] if arm=='rank_only' else each.mean()

def train(outer):
 m.check();folder=m.A/'models'/outer;folder.mkdir(parents=True,exist_ok=True)
 if (folder/'complete.json').exists():return
 tab=m.training_table(outer);tr=~tab['val'];va=tab['val'];assert tr.sum()>100 and va.sum()>100
 med,scale=m.normalizer(tab['x'][tr]);x=m.transform(tab['x'],med,scale);tab['shuffle_available']=np.ones(len(x),bool)
 shy,perm=m.shuffled_labels(tab);device='cuda' if torch.cuda.is_available() else 'cpu'
 m.write(folder/'training_table.json',dict(outer=outer,pairs=len(x),training_pairs=int(tr.sum()),validation_pairs=int(va.sum()),groups=tab['records'],
  excluded_studies=[] if outer=='all' else [outer],plan_sha256=m.b.sha256(m.A/'plan.json'),feature_dimension=x.shape[1],
  rank_donor_exclusions='Whole outer plus recipient study; recipient controls allowed',shuffled_permutations=perm,device=device))
 np.savez_compressed(folder/'normalizer.npz',median=med,scale=scale)
 xt=torch.tensor(x,device=device);yt=torch.tensor(tab['y'],device=device);pt=torch.tensor(tab['prior'],device=device);wt=torch.tensor(tab['w'],device=device)
 # A simple pair of linear residual heads is an explicit architecture baseline.
 if not (folder/'linear.joblib').exists():
  fits=[]
  for k in range(2):fits.append(Ridge(alpha=100.).fit(x[tr],tab['y'][tr,k]-tab['prior'][tr,k],sample_weight=tab['w'][tr,k]))
  joblib.dump(dict(heads=fits,median=med,scale=scale),folder/'linear.joblib',compress=3)
 validation_ids=np.flatnonzero(va);# Retain every hash0 pair; never truncate later contexts from validation.
 vi=torch.tensor(validation_ids,device=device)
 for arm in (('joint',) if outer=='all' else m.ARMS):
  dest=folder/(arm+'.pt')
  if dest.exists():assert m.b.sha256(dest)==m.read(dest.with_suffix('.json'))['sha256'];continue
  torch.manual_seed(m.SEED);torch.cuda.manual_seed_all(m.SEED);rng=np.random.default_rng(m.SEED)
  net=m.Network(x.shape[1]).to(device);opt=torch.optim.AdamW(net.parameters(),lr=m.CONFIG['learning_rate'],weight_decay=m.CONFIG['weight_decay'])
  use=tr & (tab['shuffle_available'] if arm=='shuffled_joint' else True);ids=np.flatnonzero(use)
  target=torch.tensor(shy,device=device) if arm=='shuffled_joint' else yt
  best=float('inf');weights=None;history=[];best_epoch=None
  for epoch in range(1,m.CONFIG['epochs']+1):
   net.train();losses=[]
   for step in range(m.CONFIG['steps_per_epoch']):
    ix=torch.tensor(rng.choice(ids,m.CONFIG['batch_size'],replace=len(ids)<m.CONFIG['batch_size']),device=device)
    pred=net(xt[ix],pt[ix],arm);loss=loss_values(pred,target[ix],pt[ix],wt[ix],arm)
    opt.zero_grad(set_to_none=True);loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),2.);opt.step();losses.append(float(loss.detach()))
   net.eval()
   with torch.no_grad():
    val=0.;count=0
    for st in range(0,len(vi),8192):
     j=vi[st:st+8192];v=loss_values(net(xt[j],pt[j],arm),yt[j],pt[j],wt[j],arm)
     val+=float(v)*len(j);count+=len(j)
    val/=count
   history.append(dict(epoch=epoch,training_loss=float(np.mean(losses)),source_validation_loss=val))
   if val<best:best=val;weights={k:v.detach().cpu().clone() for k,v in net.state_dict().items()};best_epoch=epoch
   print('V52',outer,arm,'epoch',epoch,'source_validation',round(val,6),flush=True)
  torch.save(dict(state_dict=weights,feature_dimension=x.shape[1],arm=arm,outer=outer,median=med,scale=scale),dest)
  m.write(dest.with_suffix('.json'),dict(sha256=m.b.sha256(dest),plan_sha256=m.b.sha256(m.A/'plan.json'),outer=outer,arm=arm,best_epoch=best_epoch,history=history,
   trained_pairs=int(use.sum()),validation_pairs=len(vi),excluded_studies=[] if outer=='all' else [outer],D3_opened=False))
  del net,target;gc.collect()
 m.write(folder/'complete.json',dict(outer=outer,plan_sha256=m.b.sha256(m.A/'plan.json'),models={p.name:m.b.sha256(p) for p in folder.glob('*.pt')},linear_sha256=m.b.sha256(folder/'linear.joblib')))
 del xt,yt,pt,wt;gc.collect();torch.cuda.empty_cache()

def load(outer,arm):
 folder=m.A/'models'/outer;p=folder/(arm+'.pt');r=m.read(p.with_suffix('.json'));assert m.b.sha256(p)==r['sha256'] and r['plan_sha256']==m.b.sha256(m.A/'plan.json')
 checkpoint=torch.load(p,map_location='cpu',weights_only=False);net=m.Network(checkpoint['feature_dimension']);net.load_state_dict(checkpoint['state_dict']);net.eval()
 return net,checkpoint

def evaluate():
 m.check();out=m.A/'predictions';out.mkdir(exist_ok=True);records=[];device='cuda' if torch.cuda.is_available() else 'cpu'
 for source in m.SOURCES:
  outer=m.GROUP[source];folder=m.A/'models'/outer;dest=out/(source+'.npz');side=dest.with_suffix('.json')
  if side.exists():
   rr=m.read(side);assert m.b.sha256(dest)==rr['sha256'];records.extend(rr['records']);continue
  models={a:load(outer,a) for a in m.ARMS}
  for net,ck in models.values():net.to(device)
  linear=joblib.load(folder/'linear.joblib');z=m.labels(source);predmean=np.zeros((300,len(z['genes'])),np.float32);predrank=predmean.copy();prior_mean=predmean.copy();prior_rank=predmean.copy()
  for ti in np.flatnonzero(z['valid'].any(1)):
   pos=np.flatnonzero(z['valid'][ti]);xx,prior=m.features(source,outer,int(ti),pos);ys=np.column_stack([z['mean'][ti].mean(0)[pos],z['rank'][ti].mean(0)[pos]])
   outputs={'v24_pooled_rank':prior.copy(),'zero_rank':np.column_stack([prior[:,0],np.zeros(len(pos))])}
   scaled=m.transform(xx,linear['median'],linear['scale']);lp=prior+np.column_stack([f.predict(scaled) for f in linear['heads']]);lp[:,0]=np.clip(lp[:,0],-3,3);lp[:,1]=np.clip(lp[:,1],-.5,.5);outputs['linear']=lp
   for arm,(net,ck) in models.items():
    xv=torch.tensor(m.transform(xx,ck['median'],ck['scale']),device=device);pv=torch.tensor(prior,device=device)
    with torch.no_grad():outputs[arm]=net(xv,pv,arm).cpu().numpy()
   predmean[ti,pos]=outputs['joint'][:,0];predrank[ti,pos]=outputs['joint'][:,1];prior_mean[ti,pos]=prior[:,0];prior_rank[ti,pos]=prior[:,1]
   validtarget=m.target_split(z['targets'][ti])==0
   for arm,pred in outputs.items():
    realrank=abs(ys[:,1])>=.02;realmean=abs(ys[:,0])>=.25
    records.append(dict(source=source,study=outer,target=str(z['targets'][ti]),target_index=int(ti),hash0_target=validtarget,arm=arm,genes=len(pos),
     mean_mse=float(np.mean((pred[:,0]-ys[:,0])**2)),rank_mse=float(np.mean((pred[:,1]-ys[:,1])**2)),
     mean_mae=float(np.mean(abs(pred[:,0]-ys[:,0]))),rank_mae=float(np.mean(abs(pred[:,1]-ys[:,1]))),
     mean_sign_accuracy=float(np.mean(np.sign(pred[realmean,0])==np.sign(ys[realmean,0]))) if realmean.any() else None,
     rank_sign_accuracy=float(np.mean(np.sign(pred[realrank,1])==np.sign(ys[realrank,1]))) if realrank.any() else None,
     mean_view_noise=float(np.mean((z['mean'][ti,0,pos]-z['mean'][ti,1,pos])**2)/4),
     rank_view_noise=float(np.mean((z['rank'][ti,0,pos]-z['rank'][ti,1,pos])**2)/4)))
  np.savez_compressed(dest,mean=predmean,rank=predrank,baseline_mean=prior_mean,baseline_rank=prior_rank,valid=z['valid'],genes=z['genes'],targets=z['targets'])
  rr=[r for r in records if r['source']==source];m.write(side,dict(source=source,excluded_study=outer,sha256=m.b.sha256(dest),plan_sha256=m.b.sha256(m.A/'plan.json'),records=rr,
   scope='All source-valid output genes per held-study target; raw mean/rank prediction only. No generated cells or official VCC scores.'))
  print('V52 evaluated',source,len(rr),'arm-target comparisons',flush=True)
  del models;gc.collect();torch.cuda.empty_cache()
 summaries={}
 for arm in ['v24_pooled_rank','zero_rank','linear',*m.ARMS]:
  study={}
  for group in m.OUTERS:
   source_stats=[]
   for source in m.SOURCES:
    if m.GROUP[source]!=group:continue
    vals=[r for r in records if r['source']==source and r['arm']==arm]
    source_stats.append({k:float(np.mean([r[k] for r in vals if r[k] is not None])) for k in ['mean_mse','rank_mse','mean_mae','rank_mae','mean_sign_accuracy','rank_sign_accuracy']})
   study[group]={k:float(np.mean([v[k] for v in source_stats])) for k in source_stats[0]}
  summaries[arm]=dict(studies=study,equal_study={k:float(np.mean([v[k] for v in study.values()])) for k in study['H1']})
 base=summaries['v24_pooled_rank']['equal_study'];joint=summaries['joint']['equal_study'];checks=dict(mean_not_worse=joint['mean_mse']<=base['mean_mse'],rank_gain_5pct=joint['rank_mse']<=.95*base['rank_mse'],
  beats_zero_rank=joint['rank_mse']<summaries['zero_rank']['equal_study']['rank_mse'],
  beats_shuffled_both=all(joint[k]<summaries['shuffled_joint']['equal_study'][k] for k in ['mean_mse','rank_mse']))
 m.write(m.A/'evaluation.json',dict(summary=summaries,checks=checks,pilot_success=all(checks.values()),scope='Repeatedly inspected public studies, newly extracted rank labels. Whole-study exclusions, not a blind private score or count-emission result.',records=len(records),D3_opened=False))
 print(json.dumps({'checks':checks,'summary':{a:v['equal_study'] for a,v in summaries.items()}},indent=2),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','train','evaluate']);p.add_argument('--source',choices=m.SOURCES);p.add_argument('--outer',choices=[*m.OUTERS,'all']);p.add_argument('--workers',type=int,default=1);args=p.parse_args()
 torch.set_num_threads(2);torch.use_deterministic_algorithms(True)
 with threadpool_limits(limits=2):
  if args.stage=='freeze':m.freeze()
  elif args.stage=='prepare':
   if args.source:m.prepare(args.source)
   elif args.workers>1:
    with ProcessPoolExecutor(max_workers=args.workers) as pool:list(pool.map(prepare_worker,m.SOURCES))
   else:
    for source in m.SOURCES:m.prepare(source)
  elif args.stage=='train':
   for outer in ([args.outer] if args.outer else [*m.OUTERS,'all']):train(outer)
  else:evaluate()
def prepare_worker(source):
 torch.set_num_threads(2)
 with threadpool_limits(limits=2):m.prepare(source)
if __name__=='__main__':main()
