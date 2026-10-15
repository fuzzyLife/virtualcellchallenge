"""Independent input, fitted-loss and prediction audits for the frozen V30 pilot."""
from pathlib import Path
import json,hashlib
import numpy as np
import h5py
from threadpoolctl import threadpool_limits
import v30_methods as m
import build_submission as b

def original_rows_audit(name):
 rel='reference/K562_gwps_raw_bulk_01.h5ad' if name=='K562_full' else 'reference/v13_resources/rpe1_raw_bulk_01.h5ad';z=np.load(m.A/'sources'/f'{name}.npz',allow_pickle=False);genes=m.read_cache()['genes'];names=z['targets'];positions=z['positions'];multi=np.flatnonzero(z['constructs']>1);selected=sorted(set([0,len(names)//2,len(names)-1]+multi[:2].tolist()));worst=0.
 with h5py.File(m.ROOT/rel,'r') as f:
  labels=f['obs/gene_transcript'].asstr()[:];n=f['obs/num_cells_filtered'][:];fold=f['obs/fold_expr'][:];core=f['obs/core_control'][:];v=f['var/gene_name'];syms=f[v.attrs['categories']].asstr()[:][v[:]];target=np.array([s.split('_',2)[1] if s.split('_',1)[0].isdigit() and len(s.split('_',2))==3 else '' for s in labels]);cols=np.flatnonzero(np.isin(syms,genes));controls=np.flatnonzero(core&(target=='non-targeting')&(n>=30)&np.isfinite(n));control=np.average(f['X'][controls,:][:,cols],axis=0,weights=n[controls]);control_total=control.sum();mapcols=[np.flatnonzero(syms[cols]==g) for g in genes[positions]];cp=np.array([control[j].sum() for j in mapcols])*1e6/control_total
  np.testing.assert_allclose(cp,z['profile'][positions],rtol=1e-12,atol=1e-10)
  for i in selected:
   ids=np.flatnonzero((target==names[i])&(n>=30)&np.isfinite(n)&np.isfinite(fold)&(fold<=.5));mean=np.average(f['X'][ids,:][:,cols],axis=0,weights=n[ids]);cpm=np.array([mean[j].sum() for j in mapcols])*1e6/mean.sum();expected=np.clip(np.log2((cpm+.1)/(cp+.1)),-3,3).astype(np.float32);error=float(abs(expected-z['y'][i]).max());assert error<1e-6;worst=max(worst,error)
 return {'source':name,'original_targets_reconstructed':len(selected),'maximum_LFC_error':worst,'normalization':'All shared competition-measured source genes before the baseline>5 output filter; not rounded pseudobulk and not whole-transcriptome single-cell counts.'}

def model_audit(name,variant):
 file=m.A/'models'/f'{name}_{variant}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(m.A/'plan.json');s=np.load(m.A/'sources'/f'{name}.npz',allow_pickle=False);model=np.load(file,allow_pickle=False)
 if rec.get('abstained'):assert not model['present'].any();return {'source':name,'variant':variant,'abstained':True}
 keys,e=m.embeddings();lookup={t:i for i,t in enumerate(keys)}
 if variant=='random_full':e=e[np.random.default_rng(m.SEED).permutation(len(e))]
 train=rec['training_targets'];rowmap={t:i for i,t in enumerate(s['targets'])};rows=np.array([rowmap[t] for t in train]);x=(e[[lookup[t] for t in train]]-model['feature_mean'])/model['feature_scale'];y=s['construct_y' if variant=='construct_full' else 'y'][rows].astype(float);w=s['n'][rows].astype(float);w=w/(w+100);w/=w.mean();gp={str(g):i for i,g in enumerate(m.read_cache()['genes'][s['positions']])};own=np.array([gp.get(t,-1) for t in train]);ok=np.flatnonzero(own>=0);mask=np.ones(y.shape,bool);mask[ok,own[ok]]=False;y[~mask]=0;bias=(w[:,None]*y).sum(0)/(w[:,None]*mask).sum(0);np.testing.assert_allclose(bias,model['bias'],atol=1e-12);center=(y-bias)/model['scale'];center[~mask]=0
 g=model['basis'][:,:min(16,model['basis'].shape[1])];ww=model['weight_r16_l1'];grad=np.zeros_like(ww);signal=np.zeros_like(ww)
 for lo in range(0,len(x),256):
  end=min(lo+256,len(x));res=x[lo:end]@ww@g.T-center[lo:end];res*=mask[lo:end];grad+=x[lo:end].T@((w[lo:end,None]*res)@g);signal+=x[lo:end].T@((w[lo:end,None]*center[lo:end])@g)
 grad+=len(x)/64*ww;relative=float(np.linalg.norm(grad)/max(np.linalg.norm(signal),1e-12));assert relative<2e-6,(name,variant,relative)
 present=model['present'];query=m.read_cache()['targets'];qx=(e[[lookup[t] for t in query[present]]]-model['feature_mean'])/model['feature_scale'];pred=np.clip((qx@ww@g.T)*model['scale']+model['bias'],-3,3)
 for ti,t in enumerate(query[present]):
  if t in gp:pred[ti,gp[t]]=0
 maximum=float(abs(pred.astype(np.float32)-model['r16_l1'][present]).max());assert maximum<2e-6
 return {'source':name,'variant':variant,'training_rows':len(train),'relative_masked_gradient':relative,'maximum_prediction_replay_error':maximum}

def main():
 plan=m.load_plan();sourcechecks=[]
 for path,h in plan['inputs'].items():assert b.sha256(m.ROOT/path)==h
 for name in m.SOURCES:
  p=m.A/'sources'/f'{name}.npz';rec=json.loads(p.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(p)
 sourcechecks=[original_rows_audit(n) for n in ['K562_full','RPE1_full']];models=[model_audit(n,v) for n,v in [('K562_full','string_full'),('K562_full','construct_full'),('RPE1_full','string_full'),('RPE1_full','string_panel'),('H1','random_full')]]
 screen=json.loads((m.A/'mean_screen.json').read_text());assert screen['plan_sha256']==b.sha256(m.A/'plan.json')
 for outer in m.OUTERS:
  for variant in ['string_full','construct_full']:
   inner=[r for r in screen['records'] if r['stage']=='inner' and r['outer']==outer and r['variant']==variant];groups={}
   for r in inner:assert r['test']!=outer;groups.setdefault(tuple(r['config']),[]).append(r['mse_ratio'])
   chosen=min(groups,key=lambda c:(np.mean(groups[c]),c));assert list(chosen)==screen['choices'][outer][variant]
 from unittest.mock import patch
 requested=[];original=m.model_arrays
 def spy(name,variant):requested.append(name);return original(name,variant)
 for outer in m.OUTERS:
  rows=np.flatnonzero(np.isin(m.read_cache()['row_families'],['HCT116','HEK293T'] if outer=='Orion' else [outer]));ci=m.read_cache()['row_contexts'][rows[0]];requested.clear()
  with patch.object(m,'model_arrays',spy):m.predict_graph(m.read_cache()['profiles'][ci],[outer],'string_full',16,1.)
  assert all(m.FAMILY[n]!=outer for n in requested)
 report={'passed':True,'plan_sha256':b.sha256(m.A/'plan.json'),'source_rows':sourcechecks,'fitted_loss_and_predictions':models,'nested_selection_recomputed':True,'actual_source_requests_respect_family_exclusion':True,'auditor_sha256':b.sha256(Path(__file__)),'scope':'Audits numerical/input integrity, not a claim of leaderboard performance. Whole-source normalization and unmeasured masks explicitly checked.'};m.write(m.A/'independent_audit.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
