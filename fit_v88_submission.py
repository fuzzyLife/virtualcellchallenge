"""Freeze and fit the sole V88 release candidate, only after its prospective release decision."""
from pathlib import Path
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import read,write,sha
A=m.A;ROOT=m.ROOT

def freeze():
 m.check();decision=read(A/'release_decision.json');assert decision['experimental_release_allowed'];assert not (A/'production_plan.json').exists()
 names=['fit_v88_submission.py','build_v88_submission.py','package_v88_submission.py','validate_v88_submission.py','v88_mixture_quantiles.py','v16_methods.py','build_submission.py'];write(A/'production_plan.json',dict(candidate='exact125',dependencies={n:sha(ROOT/n) for n in names},decision_sha256=sha(A/'release_decision.json'),scientific_plan_sha256=sha(A/'plan.json'),template_seed=20261011,emitter='residual10',automatic_upload=False,status='experimental_public_tradeoffs' if not decision['strong_support'] else 'public_supported_candidate'))

def fit():
 m.check();plan=read(A/'production_plan.json');assert sha(A/'release_decision.json')==plan['decision_sha256']
 for n,h in plan['dependencies'].items():assert sha(ROOT/n)==h,n
 from build_v66_submission import load
 manifest,genes,targets,parent,base=load();effects=[];inputs=[];trace={};dest=A/'production_inputs';dest.mkdir(exist_ok=True);assert not (A/'reference_model.json').exists()
 for ci,c in enumerate(manifest['contexts']):
  src=ROOT/'experiments/v66_assay_strength/analysis/production_inputs'/(c+'_query.npz')
  with np.load(src,allow_pickle=False) as z:query=z['query'].copy();bulk=z['bulk'].copy();library=z['library'].copy()
  pred,tr=m.predict(query,base[ci],bulk,'production',controls=False);effects.append(pred['exact125']);trace[c]=tr['exact125'];p=dest/(c+'_query.npz');np.savez_compressed(p,query=query,bulk=bulk,library=library,parent=base[ci],genes=np.asarray(genes,str));inputs.append(dict(context=c,file='context_'+c+'.h5ad',sha256=sha(ROOT/('context_'+c+'.h5ad')),query_sha256=sha(p),parent_query_sha256=sha(src)));print('V88 production response',c,flush=True);gc.collect()
 p=A/'reference_model_arrays.npz';np.savez_compressed(p,genes=np.asarray(genes,str),targets=np.asarray(targets,str),effects=np.asarray(effects,np.float32),parent_effects=base)
 write(A/'reference_model.json',dict(candidate='exact125',status=plan['status'],available_public_guards_passed=read(A/'release_decision.json')['strong_support'],confirmed_improvement=False,emitter='residual10',model_arrays_sha256=sha(p),scientific_plan_sha256=sha(A/'plan.json'),production_plan_sha256=sha(A/'production_plan.json'),decision_sha256=sha(A/'release_decision.json'),dependencies=plan['dependencies'],parent_V66_arrays_sha256=sha(ROOT/'experiments/v66_assay_strength/analysis/reference_model_arrays.npz'),context_original_inputs=inputs,guard_traces=trace,uploaded=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit']);a=p.parse_args()
 with threadpool_limits(limits=2):{'freeze':freeze,'fit':fit}[a.stage]()
