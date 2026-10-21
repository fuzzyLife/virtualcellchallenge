"""V87: independent BOLERO-inspired conditional low-rank response operator, without DNA weights."""
from pathlib import Path
from functools import lru_cache
import copy,gc
import numpy as np
import torch
from torch import nn
from scipy.stats import rankdata
from sklearn.utils.extmath import randomized_svd
from v57_pie import sha,read,write
from v62_conserved_response import data,invert_assay
from v15_ensemble_methods import cube,profiles
from v7_response_methods import FAMILIES
from v63_identity_guard import guard
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v87_conditional_adapter/analysis'
GROUP={'K562':'Replogle','CD4':'CD4','H1':'H1','HCT116':'Orion','HEK293T':'Orion'}
ARMS=('conditional25','constant25','conditional25_shuffled')
RANK=32;EPOCHS=100

def check():
 p=read(A/'plan.json')
 for k,h in p['code_sha256'].items():assert sha(ROOT/k)==h,k
 assert sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==p['source_cache_sha256']
 return p

def make_bank(d,outer,shuffle=False):
 """Select complete measured source outputs; exclude the recipient study before any fitting."""
 keep=np.array([i for i,f in enumerate(FAMILIES) if GROUP[str(f)]!=outer]);excluded=[str(f) for f in FAMILIES if GROUP[str(f)]==outer]
 y,v,ok,own=cube(d,excluded);y=y[keep];ok=ok[keep];pp=profiles(d)[keep]
 permutations=[]
 if shuffle:
  for i in range(len(y)):
   ix=np.roll(np.arange(len(d['targets'])),37+11*i);permutations.append(ix);y[i]=y[i,ix];ok[i]=ok[i,ix]
 panel=np.isin(d['genes'],d['targets']);axis=~panel & np.all(np.isfinite(pp)&(pp>5),axis=0)
 # A present target is a measured response, not an invented all-zero row.
 present=ok.sum(2)>=100
 for i in range(len(y)):
  assert present[i].any();axis &= ok[i,present[i]].all(0)
 axis=np.flatnonzero(axis);assert len(axis)>=512,(outer,len(axis))
 amount=.05*pp[:,None,axis];yy=np.log1p(amount*np.exp2(y[:,:,axis]))-np.log1p(amount)
 assert ok[:,:,axis][present].all() and np.isfinite(yy[present]).all()
 return dict(y=yy.astype(np.float32),present=present,axis=axis,profiles=pp[:,axis],source_indices=keep,source_names=FAMILIES[keep].astype(str),permutations=np.array(permutations,dtype=np.int64))

def context_rank(p,projection):
 """Map measured expression ranks to a fixed random sketch; no RNA pseudocounts."""
 p=np.asarray(p);one=p.ndim==1;p=np.atleast_2d(p);out=[]
 for row in p:
  good=np.isfinite(row);assert good.mean()>=.9
  f=np.zeros(len(row));r=rankdata(row[good],method='average');r-=r.mean();den=np.linalg.norm(r);assert den>0;f[good]=r/den;out.append(f@projection)
 arr=np.asarray(out,dtype=np.float32);return arr[0] if one else arr

def conditional_features(source,query):
 s=np.asarray(source);q=np.asarray(query);return np.concatenate([s,q,q-s],axis=-1).astype(np.float32)

class AdapterLinear(nn.Module):
 """Frozen linear weight plus independently generated rank-four factors."""
 def __init__(self,base,condition_dim=48,rank=4):
  super().__init__();self.base=copy.deepcopy(base)
  for p in self.base.parameters():p.requires_grad_(False)
  self.rank=rank;self.a=nn.Sequential(nn.Linear(condition_dim,32),nn.GELU(),nn.Linear(32,rank*base.in_features));self.b=nn.Sequential(nn.Linear(condition_dim,32),nn.GELU(),nn.Linear(32,base.out_features*rank));nn.init.zeros_(self.b[-1].weight);nn.init.zeros_(self.b[-1].bias)
 def factors(self,c):return self.a(c).reshape(-1,self.rank,self.base.in_features),self.b(c).reshape(-1,self.base.out_features,self.rank)
 def forward(self,x,c):
  a,b=self.factors(c);return self.base(x)+.1*torch.einsum('bor,bri,bi->bo',b,a,x)

class Backbone(nn.Module):
 def __init__(self):
  super().__init__();self.first=nn.Linear(RANK,64);self.last=nn.Linear(64,RANK);nn.init.zeros_(self.last.weight);nn.init.zeros_(self.last.bias)
 def forward(self,x,c=None):return x+self.last(torch.nn.functional.gelu(self.first(x)))

class ConditionalOperator(nn.Module):
 def __init__(self,base):
  super().__init__();self.first=AdapterLinear(base.first);self.last=AdapterLinear(base.last)
 def forward(self,x,c):return x+self.last(torch.nn.functional.gelu(self.first(x,c)),c)

def training_arrays(b,basis,scale,features):
 z=(b['y']@basis.T)/scale;xx=[];yy=[];cc=[];tt=[];ww=[];names=b['source_names'];groups=np.array([GROUP[n] for n in names])
 for i in range(len(names)):
  for j in range(len(names)):
   if groups[i]==groups[j]:continue
   ids=np.flatnonzero(b['present'][i]&b['present'][j]);assert len(ids)>0
   xx.append(z[i,ids]);yy.append(z[j,ids]);cc.append(np.tile(conditional_features(features[i],features[j]),(len(ids),1)));tt.append(ids);ww.append(np.full(len(ids),1/(len(ids)*np.sum(groups==groups[i])*np.sum(groups==groups[j]))))
 return [np.concatenate(a).astype(np.float32) for a in [xx,yy,cc]],np.concatenate(tt),np.concatenate(ww).astype(np.float32)

def optimize(net,arrays,targets,weights,validation,constant=False):
 x,y,c=[torch.from_numpy(v) for v in arrays];c=torch.zeros_like(c) if constant else c;w=torch.from_numpy(weights);fit=torch.from_numpy(~validation[targets]);valid=~fit;assert fit.any() and valid.any();params=[p for p in net.parameters() if p.requires_grad];opt=torch.optim.AdamW(params,lr=.001,weight_decay=.001);history=[];best=None;best_loss=np.inf;epoch_selected=None
 def loss(mask):return ((((net(x[mask],c[mask])-y[mask])**2).mean(1)*w[mask]).sum()/w[mask].sum())
 for epoch in range(EPOCHS+1):
  if epoch:
   net.train();opt.zero_grad();value=loss(fit);value.backward();nn.utils.clip_grad_norm_(params,1.);opt.step()
  net.eval()
  with torch.no_grad():train=float(loss(fit));val=float(loss(valid))
  history.append(dict(epoch=epoch,fit_loss=train,validation_loss=val))
  if val<best_loss:best_loss=val;best=copy.deepcopy(net.state_dict());epoch_selected=epoch
 net.load_state_dict(best);return dict(history=history,selected_epoch=epoch_selected,selected_validation_loss=best_loss,fit_pairs=int(fit.sum()),validation_pairs=int(valid.sum()),constant_condition=constant)

def fit(outer,shuffle=False):
 check();d=data();b=make_bank(d,outer,shuffle);validation=np.zeros(300,bool);validation[np.random.default_rng(87087).permutation(300)[:60]]=True
 train=b['y'][:,~validation][b['present'][:,~validation]];assert len(train)>=100
 _,sing,basis=randomized_svd(train,n_components=RANK,n_iter=5,random_state=87087);scale=np.maximum(np.sqrt(np.mean((train@basis.T)**2,axis=0)),.1).astype(np.float32)
 projection=np.random.default_rng(87087).normal(size=(len(b['axis']),16)).astype(np.float32);features=context_rank(b['profiles'],projection);arrays,ids,weights=training_arrays(b,basis,scale,features)
 torch.manual_seed(87087);base=Backbone();base_trace=optimize(base,arrays,ids,weights,validation)
 before={k:v.clone() for k,v in base.state_dict().items()};states={};traces={}
 for name,constant in [('conditional',False),('constant',True)]:
  if shuffle and constant:continue
  torch.manual_seed(87088);net=ConditionalOperator(base);trace=optimize(net,arrays,ids,weights,validation,constant)
  for key,old in [('first.base.weight',before['first.weight']),('first.base.bias',before['first.bias']),('last.base.weight',before['last.weight']),('last.base.bias',before['last.bias'])]:assert torch.equal(net.state_dict()[key],old)
  states[name]=net.state_dict();traces[name]=trace
 folder=A/'models';folder.mkdir(exist_ok=True);stem=outer+('_shuffled' if shuffle else '');npz=folder/(stem+'.npz');pt=folder/(stem+'.pt');assert not npz.exists() and not pt.exists()
 np.savez_compressed(npz,basis=basis,scale=scale,projection=projection,features=features,axis=b['axis'],validation_targets=np.flatnonzero(validation),training_targets=np.flatnonzero(~validation),singular_values=sing,source_names=b['source_names']);torch.save(dict(base=base.state_dict(),models=states),pt)
 write(npz.with_suffix('.json'),dict(sha256=sha(npz),weights_sha256=sha(pt),outer=outer,shuffled=shuffle,source_names=b['source_names'].tolist(),complete_output_genes=len(b['axis']),source_present_targets=b['present'].sum(1).tolist(),basis_training_response_rows=len(train),basis_training_targets=np.flatnonzero(~validation).tolist(),validation_targets=np.flatnonzero(validation).tolist(),base=base_trace,adapters=traces,plan_sha256=sha(A/'plan.json'),scope='Source-only target-held validation chooses epochs; fitting targets and basis exclude the held source-target partition. All recipient studies excluded; few independent contexts and reused biology remain limitations.'))
 print('V87 fit',stem,'genes',len(b['axis']),'validation', {k:round(v['selected_validation_loss'],6) for k,v in traces.items()},flush=True);del train,arrays;gc.collect()

@lru_cache(6)
def load(outer,shuffle):
 stem=outer+('_shuffled' if shuffle else '');p=A/'models'/(stem+'.npz');meta=read(p.with_suffix('.json'));assert sha(p)==meta['sha256'];assert sha(p.with_suffix('.pt'))==meta['weights_sha256']
 with np.load(p,allow_pickle=False) as z:arr={k:z[k].copy() for k in z.files}
 weights=torch.load(p.with_suffix('.pt'),map_location='cpu',weights_only=True);base=Backbone();base.load_state_dict(weights['base']);nets={}
 for key,state in weights['models'].items():net=ConditionalOperator(base);net.load_state_dict(state);net.eval();nets[key]=net
 return arr,nets

def response(query,outer,shuffle=False):
 arr,nets=load(outer,shuffle);b=make_bank(data(),outer,shuffle);assert np.array_equal(b['axis'],arr['axis']);assert np.array_equal(b['source_names'],arr['source_names']);axis=arr['axis'];feature=context_rank(query[axis],arr['projection']);source=arr['features'];similar=np.sum(source*feature,axis=1)/np.maximum(np.linalg.norm(source,axis=1)*np.linalg.norm(feature),1e-8);groups=np.array([GROUP[n] for n in b['source_names']]);w=np.exp(2*(similar-similar.max()))/np.array([np.sum(groups==g) for g in groups]);outputs={}
 for name,net in nets.items():
  total=np.zeros((300,len(axis)));den=np.zeros(300)
  for i in range(len(source)):
   ids=np.flatnonzero(b['present'][i]);z=(b['y'][i,ids]@arr['basis'].T)/arr['scale'];c=np.tile(conditional_features(source[i],feature),(len(ids),1));c*=0 if name=='constant' else 1
   with torch.no_grad():pred=net(torch.from_numpy(z),torch.from_numpy(c)).numpy()*arr['scale']
   total[ids]+=w[i]*(pred@arr['basis']);den[ids]+=w[i]
  outputs[name]=np.divide(total,den[:,None],out=np.zeros_like(total),where=den[:,None]>0)
 return outputs,den>0,axis,dict(source_names=b['source_names'].tolist(),source_weights=w.tolist(),query_measured_fraction=float(np.isfinite(query[axis]).mean()))

def predict(query,base,bulk,outer):
 check();d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);out={};trace={}
 for shuffle in [False,True]:
  arr,_=load(outer,shuffle)
  if np.isfinite(query[arr['axis']]).mean()<.9:
   for name in (['conditional25_shuffled'] if shuffle else ['conditional25','constant25']):out[name]=base.copy();trace[name]={'abstained_missing_query':True}
   continue
  responses,support,axis,tr=response(query,outer,shuffle)
  for name,yy in responses.items():
   arm=name+'25'+('_shuffled' if shuffle else '');full=np.zeros_like(base,dtype=float);full[:,axis]=yy;desired=invert_assay(full,query);mask=np.zeros_like(base,bool);mask[:,axis]=True;mask&=support[:,None]&np.isfinite(query)[None,:]&(query[None,:]>0);proposal=np.where(mask,base+.25*np.clip(desired-base,-.5,.5),base).astype(np.float32);proposal[np.arange(300),own]=0;value,diag=guard(base,proposal,bulk,own,panel);inactive=~np.any(proposal!=base,axis=1);value[inactive]=base[inactive];out[arm]=value;trace[arm]=dict(diag,routing=tr,changed_targets=int(np.any(value!=base,axis=1).sum()),corrected_axis_genes=len(axis))
 return out,trace
