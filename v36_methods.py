"""V36: two independent population views, mean and dispersion supervision.

Uses moment summaries of control sets, not observed recipient treatment means.
Raw-cell training excludes the same complete families as retrieval training.
"""
import copy,gc,json
import numpy as np
import torch
from torch import nn
from scipy import sparse
import v35_methods as m
A=m.P36

def view_stats(x,ids,axis):
 xx=x[ids].astype(float);libs=xx.sum(1).A1
 y=xx[:,axis].toarray()*1e6/libs[:,None];mu=y.mean(0);var=np.maximum(y.var(0),0)
 return mu,var

def build_populations(g,a,records):
 groups=[];receipts=[];data=m.old.old.read_cache();axis=g['axis']
 for source,rawname in m.RAW_NAME.items():
  if source not in g['sources']:continue
  selected=[i for i,r in enumerate(records) if r['source']==source]
  if not selected:continue
  x,meta,measured,rec=m.raw(rawname,data['genes']);assert measured[axis].all()
  labels=meta.gene_target.to_numpy(str);lanes=meta['sample'].to_numpy(str);ntc=np.flatnonzero(labels=='Non-Targeting');stats=m.control_stats(x[ntc]);rng=np.random.default_rng(m.SEED+m.CI[source])
  pools={s:np.flatnonzero((labels=='Non-Targeting')&(lanes==s)) for s in sorted(set(lanes[ntc]))}
  used=[]
  for ri in selected:
   t=records[ri]['target'];ids=np.flatnonzero(labels==t);ids=ids[np.isin(lanes[ids],list(pools))]
   n=min(32,len(ids)//2)
   if n<16:continue
   yviews=[];vviews=[];valids=[];traces=[]
   for vi in range(8):
    order=rng.permutation(ids)[:2*n];left,right=order[:n],order[n:]
    ctrls=[[],[]]
    for lane in sorted(set(lanes[order])):
     pool=pools[lane];nc=min(64,max(8,int(np.sum(lanes[order]==lane))),len(pool)//2)
     if nc<4:continue
     cc=rng.permutation(pool)[:2*nc];ctrls[0].extend(cc[:nc].tolist());ctrls[1].extend(cc[nc:].tolist())
    assert ctrls[0] and ctrls[1] and not set(ctrls[0])&set(ctrls[1]) and not set(left)&set(right)
    yy=[];vv=[];mm=[]
    for take,ctrl in zip([left,right],ctrls):
     # Match control lane proportions to this particular treated view.
     tm,tv=view_stats(x,take,axis);cm=np.zeros(len(axis));cv=cm.copy()
     for lane in set(lanes[take]):
      cc=np.array([i for i in ctrl if lanes[i]==lane]);assert len(cc)
      ll,lv=view_stats(x,cc,axis);w=float(np.mean(lanes[take]==lane));cm+=w*ll;cv+=w*lv
     yy.append(np.clip(np.log2((tm+.1)/(cm+.1)),-3,3))
     ratio=np.log2((tv/(tm+1)**2+.01)/(cv/(cm+1)**2+.01));vv.append(np.clip(ratio,-1.5,1.5))
     mm.append((cm>5)&(data['genes'][axis]!=t))
    yviews.append(yy);vviews.append(vv);valids.append(mm)
    trace={'view':vi,'treated_sizes':[len(left),len(right)],'control_sizes':[len(v) for v in ctrls],'disjoint':True}
    if vi==0:trace.update(treated_rows=[left.tolist(),right.tolist()],control_rows=ctrls)
    traces.append(trace)
   pp,var,drop=stats
   groups.append({'row':ri,'y':np.array(yviews,np.float32),'dispersion':np.array(vviews,np.float32),'valid':np.array(valids,bool),'control_var':np.clip(np.log1p(var[axis]/np.maximum(pp[axis]**2,1e-14))/5,0,2).astype(np.float32),'dropout':drop[axis]})
   used.append({'training_row':ri,'target':t,'validation':records[ri]['validation'],'views':traces})
  receipts.append({'source':source,'family':m.family(source),'raw_source':rawname,'source_record':rec,'control_rows':ntc.tolist(),'groups':used})
  del x;gc.collect();print('V36 populations',g['outer'],source,len(used),flush=True)
 assert groups
 return groups,receipts

def population_batch(a,pop,records,ids,cols,rng,device):
 rows=[pop[i]['row'] for i in ids];b,ge=m.tensor_batch(a,rows,cols,device);ys=[];vs=[];valid=[]
 for j,i in enumerate(ids):
  p=pop[i];v=int(rng.integers(8));ys.append(p['y'][v,:,cols].T);vs.append(p['dispersion'][v,:,cols].T);valid.append(p['valid'][v,:,cols].T)
  b['control_var'][j]=torch.tensor(p['control_var'][cols],device=device);b['dropout'][j]=torch.tensor(p['dropout'][cols],device=device);b['moments_known'][j]=1.
 # Advanced indexing produces (genes,2); transpose above gives (2,genes).
 return b,ge,torch.tensor(np.array(ys),device=device),torch.tensor(np.array(vs),device=device),torch.tensor(np.array(valid),device=device,dtype=torch.float32)

def population_loss(net,b,ge,y,variance,valid):
 pred,v,delta=net(b,ge);mask=valid.prod(1);ymean=y.mean(1);noise=(y[:,0]-y[:,1])**2
 reliability=torch.clamp(1/(.1+noise),min=.1,max=5)*mask
 reliability=reliability/torch.clamp(reliability.mean(1,keepdim=True),min=1e-6)
 # Both independent population views supervise one mean, rather than pairing cells.
 den=torch.clamp((((b['base'][:,None,:]-y)**2)*reliability[:,None,:]).sum((1,2)),min=.04*mask.sum(1))
 mean=((((pred[:,None,:]-y)**2)*reliability[:,None,:]).sum((1,2))/den).mean()
 vl=nn.functional.smooth_l1_loss(v,variance.mean(1),reduction='none');vl=(vl*mask).sum()/torch.clamp(mask.sum(),min=1)
 return mean+.1*vl+.05*(delta*delta).mean(),mean,vl

def validate_pop(net,a,pop,records,device):
 groups=[i for i,p in enumerate(pop) if records[p['row']]['validation']];rng=np.random.default_rng(m.SEED);families={}
 with torch.no_grad():
  for f in sorted({records[pop[i]['row']]['family'] for i in groups}):
   ids=[i for i in groups if records[pop[i]['row']]['family']==f];vals=[]
   for st in range(0,len(ids),8):
    cols=np.arange(0,len(a['gene']),max(1,len(a['gene'])//1024));b,ge,y,v,mask=population_batch(a,pop,records,ids[st:st+8],cols,rng,device);loss,mean,var=population_loss(net,b,ge,y,v,mask);vals.append(float(loss))
   families[f]=float(np.mean(vals))
 return float(np.mean(list(families.values()))),families

def train(outer):
 m.setup();m.check();torch.manual_seed(m.SEED+1);torch.cuda.manual_seed_all(m.SEED+1);rng=np.random.default_rng(m.SEED+1)
 net,g,parent=m.load(outer,35);a,rec=m.table(g);pop,receipt=build_populations(g,a,rec)
 folder=A/'populations';folder.mkdir(parents=True,exist_ok=True);m.write(folder/(outer+'.json'),{'plan_sha256':m.sha(m.A/'plan.json'),'outer':outer,'sources':receipt})
 tr=np.array([i for i,r in enumerate(rec) if not r['validation']]);ptr=np.array([i for i,p in enumerate(pop) if not rec[p['row']]['validation']]);pr=[rec[p['row']] for p in pop]
 device='cuda' if torch.cuda.is_available() else 'cpu';net.to(device);opt=torch.optim.AdamW(net.parameters(),lr=.0002,weight_decay=.001);history=[];score=float('inf');best=None;weights=None
 for epoch in range(1,m.EPOCHS+1):
  net.train();losses=[]
  for step in range(m.STEPS):
   cols=rng.choice(len(g['axis']),min(512,len(g['axis'])),replace=False);ids=m.balanced_rows(rec,tr,8,rng);b,ge=m.tensor_batch(a,ids,cols,device);p,_,delta=net(b,ge);summary=m.loss_mean(p,b['truth'],b['base'],b['valid'])
   ids=m.balanced_rows(pr,ptr,8,rng);b,ge,y,v,mask=population_batch(a,pop,rec,ids,cols,rng,device);pl,_,_=population_loss(net,b,ge,y,v,mask);loss=.5*summary+.5*pl+.025*(delta*delta).mean();opt.zero_grad(set_to_none=True);loss.backward();nn.utils.clip_grad_norm_(net.parameters(),1.);opt.step();losses.append(float(loss.detach()))
  net.eval();summary,fams=m.validate(net,a,rec,device);population,pf=validate_pop(net,a,pop,rec,device);val=.5*summary+.5*population;history.append({'epoch':epoch,'loss':float(np.mean(losses)),'validation_MSE_ratio':summary,'validation_population_loss':population,'selection_loss':val,'families':fams,'population_families':pf})
  if val<score:score=val;best=epoch;weights=copy.deepcopy(net.state_dict())
  print('V36',outer,'epoch',epoch,'validation',val,'summary',summary,flush=True)
 net.load_state_dict(weights);m.save(net,g,rec,history,best,36,{'parent_model_sha256':parent['sha256'],'population_manifest_sha256':m.sha(folder/(outer+'.json')),'population_groups':len(pop),'population_fit_targets':sorted({pr[i]['target'] for i in ptr}),'population_validation_targets':sorted({r['target'] for r in pr if r['validation']})});del a,pop,net;gc.collect();torch.cuda.empty_cache()


def emit(template,effect,own,seed,prepared,atlas,dispersion=None):
 from v25_methods import emit as old_emit
 if dispersion is None:return old_emit(template,effect,own,seed,prepared,atlas)
 libs,p,bulk=prepared;libs=np.asarray(libs,np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();expected=libs[:,None]*wp;emp=template.toarray().astype(float)*factor;emp*=libs[:,None]/emp.sum(1,keepdims=True)
 alpha=np.clip(.1*np.exp2(np.clip(dispersion,-1.5,1.5)/2),.05,.25)
 initial=(.35-alpha[None,:])*expected+.65*np.random.default_rng(seed).poisson(expected)+alpha[None,:]*emp;initial/=initial.sum(1,keepdims=True)
 count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200);np.testing.assert_array_equal(count.sum(1,dtype=np.int64),libs)
 assert count.dtype.kind in 'iu' and count.min()>=0
 rp=(count/libs[:,None]).mean(0);rb=count.sum(0,dtype=np.int64)/libs.sum()
 return sparse.csr_matrix(count),{'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'empirical_weight_min':float(alpha.min()),'empirical_weight_max':float(alpha.max())}