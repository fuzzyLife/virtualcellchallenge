"""Diagnose V26 source-control posterior/prior moments without changing any predictor."""
from pathlib import Path
import json,argparse
import numpy as np,torch
from threadpoolctl import threadpool_limits
from v26_network import Decoder
from v26_methods import ROOT,A,KINDS,load_plan,load_decoder,write
import build_submission as b

def main(key):
 load_plan();file=A/'training'/f'{key}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file)
 with np.load(file,allow_pickle=False) as z:d={n:z[n].copy() for n in z.files}
 labels=np.array([g['label'] for g in rec['group_records']]);ids=np.flatnonzero((~d['fit'])&(labels[d['groups']]=='Non-Targeting'));assert len(ids)>20;out={}
 for kind in KINDS:
  model=load_decoder(key,kind);net=Decoder(d['features'].shape[1],18533);net.load_state_dict({n:torch.from_numpy(v) for n,v in model['weights'].items()});net.eval();mu=[];lv=[];states=d['state'] if kind!='state_shuffled' else d['state'][d['shuffled']];states=np.zeros_like(states) if kind=='independent' else states
  with torch.no_grad():
   for lo in range(0,len(ids),256):
    take=ids[lo:lo+256];m,l=net.posterior(torch.from_numpy(d['features'][take]),torch.from_numpy(states[take]),torch.from_numpy(d['conditions'][d['groups'][take]]));mu.append(m.numpy());lv.append(l.numpy())
  mu=np.concatenate(mu);lv=np.concatenate(lv);mean=mu.mean(0);variance=np.exp(lv).mean(0)+mu.var(0);cross=np.corrcoef(np.c_[mu,states[ids]].T)[:16,16:] if kind!='independent' else np.zeros((16,8));out[kind]={'model_sha256':model['record']['sha256'],'validation_control_cells':len(ids),'aggregated_posterior_mean':mean.tolist(),'aggregated_posterior_total_variance':variance.tolist(),'mean_norm':float(np.linalg.norm(mean)),'variance_min':float(variance.min()),'variance_median':float(np.median(variance)),'variance_max':float(variance.max()),'diagonal_Gaussian_moment_KL_proxy':float(.5*np.sum(mean*mean+variance-np.log(variance)-1)),'maximum_state_posterior_mean_correlation':float(np.nanmax(abs(cross)))}
 write(A/f'{key}_latent_diagnostic.json',{'plan_sha256':b.sha256(A/'plan.json'),'training_arrays_sha256':b.sha256(file),'models':out,'diagnostic_sha256':b.sha256(Path(__file__)),'scope':'Descriptive source-control held-cell posterior moments versus the N(0,I) generation prior. Aggregated posterior is a mixture; Gaussian diagonal proxy is not its exact KL. Correlation/mismatch does not identify a causal mechanism or establish unseen-context performance. No parameters or release checks changed.'});print(json.dumps({k:{n:v[n] for n in ['validation_control_cells','mean_norm','variance_min','variance_median','variance_max','maximum_state_posterior_mean_correlation']} for k,v in out.items()},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--family',default='H1');a=p.parse_args();torch.set_num_threads(2)
 with threadpool_limits(limits=2):main(a.family)
