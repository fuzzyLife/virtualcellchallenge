"""One-step Stack conditional mean readout; exact author mixing/forward, no NB sampling."""
import math
import numpy as np,torch

def mixed_windows(base,query,n_cells=512,prompt=.25,context=.4,seed=20020):
 n_test=max(1,int(n_cells*(1-prompt-context)));n_base=n_cells-n_test;rng=np.random.RandomState(seed);order=rng.permutation(base.shape[0]);windows=[];keep=[];ptr=0
 for i in range(math.ceil(query.shape[0]/n_test)):
  ids=np.arange(i*n_test,min((i+1)*n_test,query.shape[0]));real=len(ids)
  if real<n_test:ids=np.r_[ids,np.arange(n_test-real)]
  source=order[np.arange(ptr,ptr+n_base)%len(order)];ptr=(ptr+n_base)%len(order);windows.append(np.vstack([base[source].toarray(),query[ids].toarray()]).astype(np.float32));keep.extend([False]*n_base+[True]*real+[False]*(n_test-real))
 return np.stack(windows),np.asarray(keep,bool)

@torch.inference_mode()
def predict_means(model,base,query,seed=20020,batch_size=4):
 windows,keep=mixed_windows(base,query,model.n_cells,seed=seed);means=[];disp=[];device=next(model.parameters()).device
 for lo in range(0,len(windows),batch_size):
  features=torch.from_numpy(windows[lo:lo+batch_size]).to(device);observed=features.sum(-1,keepdim=True);tokens=model._reduce_and_tokenize(torch.log1p(features));n_kept=int(.25*model.n_cells);n_context=int((.25+.4)*model.n_cells);expanded=torch.zeros_like(tokens);expanded[:,n_context:]=1;tokens=tokens+model.query_pos_embedding[None,None]*expanded;mask=torch.zeros(model.n_cells,model.n_cells,dtype=torch.bool,device=device);mask[:n_kept,n_kept:]=True;x=model._run_attention_layers(tokens,gene_attn_mask=mask[None,None]);embedding=x.reshape(features.shape[0],model.n_cells,-1);mu,theta,_=model._compute_nb_parameters(embedding,observed);means.append(mu.cpu().numpy().reshape(-1,model.n_genes));disp.append(theta.cpu().numpy().reshape(-1,model.n_genes))
 return np.concatenate(means)[keep],np.concatenate(disp)[keep]
