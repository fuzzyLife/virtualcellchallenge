"""Fit V63 with the fixed assay predictor plus an outcome-free identity guard."""
import argparse,gc
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import v63_identity_guard as m
from v62_pool_cache import install
from v57_pie import sha,read,write
from v6_methods import moments
install();ROOT=m.ROOT;A=m.A

def freeze():
 assert not (A/'production_plan.json').exists();m.check();m.parent.check()
 files=['fit_v63_submission.py','build_v63_submission.py','package_v63_submission.py','validate_v63_submission.py','verify_v63_effects.py','v63_identity_guard.py','v62_conserved_response.py','v62_pool_cache.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py']
 write(A/'production_plan.json',dict(candidate='guarded_assay',scientific_plan_sha256=sha(A/'plan.json'),dependencies={n:sha(ROOT/n) for n in files},emitter='residual10',template_seed=20261011,automatic_upload=False,D3_opened=False,selection='Only frozen guarded primary, with source-target shuffled diagnostic. No alternative promotion after outcomes. Explicit --experimental mode records failed gates and cannot claim validated improvement.'))

def fit(experimental=False):
 m.check();m.parent.check();p=read(A/'production_plan.json');assessment=read(A/'count_assessment.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert assessment['complete'] and (experimental or assessment['release_supported'])
 from build_v24_submission import load
 import build_submission as b
 manifest,genes,targets,oldmodel,base=load();own=np.array([list(genes).index(t) for t in targets]);panel=np.isin(genes,targets);effects=[];traces=[];changes=[];pdir=A/'production_inputs';pdir.mkdir(exist_ok=True);routing={}
 assert not (A/'reference_model.json').exists()
 for ci,c in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(ROOT,c,genes,manifest);depth,mean,bulk,_=moments(x);q=mean*1e6
  candidate=m.parent.predict(q,base[ci],())['assay_center50'];e,trace=m.guard(base[ci],candidate,bulk,own,panel);effects.append(e);routing[c]=trace;path=pdir/(c+'_query.npz')
  np.savez_compressed(path,query=q,query_mean=q,baseline=base[ci],unprotected=candidate,genes=np.asarray(genes,str),bulk=bulk,library=library)
  traces.append(dict(context=c,file=f'context_{c}.h5ad',sha256=sha(ROOT/f'context_{c}.h5ad'),query_sha256=sha(path),control_cells=x.shape[0]))
  for ti,t in enumerate(targets):changes.append(dict(context=c,target_gene=str(t),sum_absolute_LFC_correction=float(abs(e[ti]-base[ci,ti]).sum()),changed_genes=int(np.count_nonzero(e[ti]!=base[ci,ti]))))
  del x;gc.collect();print('V63 production',c,'full correction',trace['full_correction_targets'],flush=True)
 path=A/'reference_model_arrays.npz';np.savez_compressed(path,genes=np.asarray(genes,str),targets=np.asarray(targets,str),effects=np.asarray(effects,np.float32),parent_effects=base)
 pd.DataFrame(changes).to_csv(A/'production_response_changes.csv',index=False)
 write(A/'reference_model.json',dict(model='V63 identity-guarded assay-space response transfer',candidate='guarded_assay',status='public_supported_candidate' if assessment['release_supported'] else 'experimental_failed_public_guards',available_public_guards_passed=assessment['release_supported'],confirmed_improvement=False,emitter='residual10',model_arrays_sha256=sha(path),scientific_plan_sha256=sha(A/'plan.json'),production_plan_sha256=sha(A/'production_plan.json'),count_assessment_sha256=sha(A/'count_assessment.json'),context_original_inputs=traces,dependencies=p['dependencies'],parent_V24_arrays_sha256=sha(ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),guard_traces=routing,D3_opened=False,uploaded=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit']);p.add_argument('--experimental',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  else:fit(a.experimental)
