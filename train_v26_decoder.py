"""Prepare bounded source-only V26 training samples and fit three prespecified decoders."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
from pathlib import Path
import argparse,json,gc,copy,time
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits
import torch
from v26_network import Decoder
from v26_methods import ROOT,A,FAMILIES,SOURCES,KINDS,load_plan,family_key,raw,geometry,descriptor,features,write
from fit_v7_response import read_cache
import build_submission as b

def prepare(key):
 plan=load_plan();folder=A/'training';folder.mkdir(exist_ok=True);dest=folder/f'{key}.npz'
 if dest.exists():
  rec=json.loads(dest.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(dest) and rec['plan_sha256']==b.sha256(A/'plan.json');return
 data=read_cache();genes=data['genes'];targets=data['targets'];geo=geometry(key);training=[f for f in SOURCES if family_key(f)!=key];rowkeep=np.isin(data['row_families'],[family_key(f) for f in training]);v=np.where(data['valid'][rowkeep],abs(data['truth'][rowkeep]),np.nan)
 with np.errstate(invalid='ignore'):sensitivity=np.nanmedian(v,axis=0)
 sensitivity=np.nan_to_num(sensitivity,nan=100);ci={'H1':4,'HCT116':5,'HEK293T':6,'CD4_Rest_D2':1,'CD4_Stim48_D2':3};common=np.isfinite(data['profiles'][[ci[f] for f in training]]).all(0)&~np.isin(genes,targets);eligible=np.flatnonzero(common&(np.nanmin(data['profiles'][[ci[f] for f in training]],axis=0)>5));state_axis=eligible[np.argsort(sensitivity[eligible],kind='stable')[:512]];assert len(state_axis)==512
 controls=[];blocks=[];records=[];group_base=[];group_control=[];group_condition=[];allgroups=[];allfit=[];allstatesraw=[];source_index=[];output_mask=np.ones(len(genes),bool);offset=0
 for fi,f in enumerate(training):
  x,m,measured,source=raw(f,genes);output_mask&=measured;labels=m.gene_target.to_numpy(str);batch=m['sample'].to_numpy(str);rng=np.random.default_rng(26000+fi);selected=[];fit=[];batchbase={};batchcenter={};control_rows=[]
  for lane in sorted(set(batch[labels=='Non-Targeting'])):
   ids=np.flatnonzero((labels=='Non-Targeting')&(batch==lane));ids=rng.permutation(ids)[:64]
   if len(ids)<10:continue
   flags=np.arange(len(ids))>=max(2,len(ids)//10);selected.extend(ids.tolist());fit.extend(flags.tolist());fitrows=ids[flags];d=np.asarray(x[fitrows].sum(1)).ravel();pp=np.asarray(x[fitrows].multiply((1/d)[:,None]).mean(0)).ravel();batchbase[lane]=pp;batchcenter[lane]=np.log1p(1e4*pp[state_axis]);control_rows.extend(fitrows.tolist())
  assert batchbase
  for ti,t in enumerate(targets):
   ids=np.flatnonzero((labels==t)&np.isin(batch,list(batchbase)));ids=rng.permutation(ids)[:96]
   if len(ids)<16:continue
   flags=np.arange(len(ids))>=max(2,len(ids)//10);selected.extend(ids.tolist());fit.extend(flags.tolist())
  selected=np.asarray(selected,int);fit=np.asarray(fit,bool);xx=x[selected].tocsr();ll=labels[selected];bb=batch[selected];libs=np.asarray(xx.sum(1)).ravel();assert (libs>0).all();norm=xx.multiply((1/libs)[:,None]).tocsr();state=np.log1p(1e4*norm[:,state_axis].toarray())-np.array([batchcenter[a] for a in bb]);controls.append(state[(ll=='Non-Targeting')&fit]);allstatesraw.append(state.astype(np.float32));gids=np.zeros(len(selected),int)
  for label in sorted(set(ll)):
   keep=(ll==label);train=keep&fit;assert train.any();pp=np.asarray(norm[train].mean(0)).ravel();ctrl=np.mean([batchbase[a] for a in bb[train]],axis=0);effect=np.clip(np.log2((pp+1e-8)/(ctrl+1e-8)),-3,3) if label!='Non-Targeting' else np.zeros(len(genes));effect[np.isin(genes,targets)]=0;gids[keep]=len(group_base);smooth=pp+.001*ctrl+1e-8*measured;smooth/=smooth.sum();group_base.append(smooth.astype(np.float32));group_control.append(ctrl.astype(np.float32));group_condition.append(descriptor(effect,ctrl,geo));allgroups.append({'context':f,'label':str(label),'training_cells':int(train.sum()),'validation_cells':int((keep&~fit).sum())})
  blocks.append(xx);allfit.append(fit);source_index.extend([fi]*len(selected));records.append({'context':f,'family':family_key(f),'source':source,'selected_original_rows':selected.tolist(),'fit_flags':fit.tolist(),'control_fit_rows':control_rows,'projected_libraries':libs.astype(int).tolist(),'original_libraries':m.total_counts.to_numpy(float)[selected].tolist(),'batch_labels':bb.tolist(),'group_ids':gids.tolist()});offset+=len(selected);del x,xx,norm;gc.collect();print('V26 training rows selected',key,f,len(selected),flush=True)
 controls=np.concatenate(controls);scale=np.maximum(controls.std(0),.2);_,_,vt=randomized_svd(controls/scale,n_components=8,random_state=26026);state_basis=vt.T;state=np.clip(np.concatenate(allstatesraw)/scale@state_basis,-4,4).astype(np.float32);matrix=sparse.vstack(blocks,format='csr');fit=np.concatenate(allfit);groups=np.concatenate([np.array(r['group_ids']) for r in records]);sources=np.asarray(source_index);shuffled=np.arange(len(fit));start=0
 for rec in records:
  batch=np.array(rec['batch_labels']);n=len(batch)
  for lane in set(batch):
   for flag in [False,True]:
    ids=np.flatnonzero((batch==lane)&(fit[start:start+n]==flag))+start;shuffled[ids]=np.random.default_rng(26026+start+int(flag)).permutation(ids)
  start+=n
 output_mask&=~np.isin(genes,targets);base=np.array(group_base);controls_base=np.array(group_control);conditions=np.array(group_condition);qfeatures=features(matrix,base[groups],geo['positions']);sparse.save_npz(folder/f'{key}_counts.npz',matrix);np.savez_compressed(dest,genes=genes.astype(str),targets=targets.astype(str),state_axis=state_axis,state_basis=state_basis,state_scale=scale,state=state,shuffled=shuffled,fit=fit,groups=groups,base=base,control=controls_base,conditions=conditions,features=qfeatures,output_mask=output_mask,sources=sources);write(dest.with_suffix('.json'),{'sha256':b.sha256(dest),'counts_sha256':b.sha256(folder/f'{key}_counts.npz'),'plan_sha256':b.sha256(A/'plan.json'),'excluded_family':key,'training_contexts':training,'rows':matrix.shape[0],'genes':matrix.shape[1],'trained_output_genes':int(output_mask.sum()),'group_records':allgroups,'original_inputs':records,'scope':'Group means use fitting cells only; posterior can encode treatment residuals. Stable-feature state is a proxy, not an identified individual pre-treatment counterfactual. Source projected-axis libraries used for multinomial likelihood; original complete libraries retained separately.'});print('V26 training data prepared',key,matrix.shape,flush=True)

def train(key,kind):
 plan=load_plan();prepare(key);file=A/'models'/f'{key}_{kind}.npz';file.parent.mkdir(exist_ok=True)
 if file.exists():
  rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file);return
 ds=A/'training'/f'{key}.npz';rec=json.loads(ds.with_suffix('.json').read_text());counts=sparse.load_npz(ds.with_name(f'{key}_counts.npz'))
 with np.load(ds,allow_pickle=False) as z:d={n:z[n].copy() for n in z.files}
 assert rec['sha256']==b.sha256(ds) and rec['counts_sha256']==b.sha256(ds.with_name(f'{key}_counts.npz'));torch.manual_seed(26026);torch.cuda.manual_seed_all(26026);torch.set_num_threads(2);torch.use_deterministic_algorithms(True);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;device='cuda' if torch.cuda.is_available() else 'cpu';net=Decoder(d['features'].shape[1],counts.shape[1]).to(device);optimizer=torch.optim.AdamW(net.parameters(),lr=.001,weight_decay=.0001);mask=torch.from_numpy(d['output_mask'].astype(np.float32)).to(device);states=d['state'] if kind!='state_shuffled' else d['state'][d['shuffled']];states=np.zeros_like(states) if kind=='independent' else states;trainids=np.flatnonzero(d['fit']);validids=np.flatnonzero(~d['fit']);rng=np.random.default_rng(26026);history=[];best=float('inf');weights=None
 def batch(ids):
  y=counts[ids].toarray().astype(np.float32);libs=y.sum(1);g=d['groups'][ids];return [torch.from_numpy(a).to(device) for a in [d['features'][ids],states[ids],d['conditions'][g],d['base'][g],y/libs[:,None]]]
 for epoch in range(12):
  net.train();losses=[];start=time.time()
  for lo in range(0,len(trainids),128):
   if lo==0:order=rng.permutation(trainids)
   x,s,c,base,y=batch(order[lo:lo+128]);optimizer.zero_grad(set_to_none=True);lp,kl=net(x,s,c,base,mask);loss=-(y*lp).sum(1).mean()+.001*kl.mean();assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),5.);optimizer.step();losses.append(float(loss.detach()))
  net.eval();vals=[]
  with torch.no_grad():
   torch.manual_seed(26126)
   for lo in range(0,len(validids),128):
    x,s,c,base,y=batch(validids[lo:lo+128]);lp,kl=net(x,s,c,base,mask);vals.append(float((-(y*lp).sum(1)+.001*kl).mean()))
  value=float(np.mean(vals));history.append({'epoch':epoch+1,'training_elbo':float(np.mean(losses)),'validation_elbo':value,'seconds':time.time()-start});print('V26 fit',key,kind,epoch+1,history[-1],flush=True)
  if value<best:best=value;weights={n:v.detach().cpu().numpy().copy() for n,v in net.state_dict().items()}
 np.savez_compressed(file,**{'w__'+n:v for n,v in weights.items()},output_mask=d['output_mask'],state_axis=d['state_axis'],state_basis=d['state_basis'],state_scale=d['state_scale']);write(file.with_suffix('.json'),{'sha256':b.sha256(file),'plan_sha256':b.sha256(A/'plan.json'),'training_record_sha256':b.sha256(ds.with_suffix('.json')),'training_arrays_sha256':b.sha256(ds),'excluded_family':key,'training_contexts':rec['training_contexts'],'kind':kind,'selection':'Minimum held-cell validation ELBO among12 fixed epochs; no treatment scorer outcomes used','history':history,'runtime':{'torch':torch.__version__,'device':device,'cuda':torch.version.cuda,'gpu':torch.cuda.get_device_name(0) if device=='cuda' else None,'deterministic_algorithms':True},'network':'4096-feature posterior ->128 ->16 latent; decoder16 latent+8 state+16 effect ->64 ->128 ->18533 bounded logit residuals; masked unmeasured/panel output residuals.'});del net,counts,d;gc.collect();torch.cuda.empty_cache();print('V26 decoder saved',key,kind,flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','train']);p.add_argument('--family',choices=FAMILIES+['all']);a=p.parse_args()
 for key in ([a.family] if a.family else FAMILIES):
  if a.stage=='prepare':prepare(key)
  else:
   for kind in KINDS:train(key,kind)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
