"""Fit only the count-tested V66 candidate, recording failed guards when experimental."""
import argparse,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import v66_assay_strength as m
from v57_pie import sha,read,write
from v6_methods import moments
ROOT=m.ROOT;A=m.A

def freeze():
 m.check();m.old.check();assert not (A/'production_plan.json').exists();selected=read(A/'mean_results.json')['selected'];assert selected
 files=['fit_v66_submission.py','build_v66_submission.py','package_v66_submission.py','validate_v66_submission.py','verify_v66_effects.py','v66_assay_strength.py','v62_conserved_response.py','v63_identity_guard.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py']
 write(A/'production_plan.json',dict(candidate=selected,scientific_plan_sha256=sha(A/'plan.json'),mean_selection_sha256=sha(A/'mean_results.json'),dependencies={n:sha(ROOT/n) for n in files},emitter='residual10',template_seed=20261011,automatic_upload=False,D3_opened=False,selection='Candidate locked solely by frozen desired-mean screen before any new count outcome; count guard failures remain explicit. No post-count promotion.'))

def fit(experimental=False):
 m.check();m.old.check();p=read(A/'production_plan.json');assessment=read(A/'count_assessment.json');assert assessment['complete'] and (assessment['release_supported'] or experimental)
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert sha(A/'mean_results.json')==p['mean_selection_sha256'];assert read(A/'boundary_audit.json')['passed']
 from build_v24_submission import load
 import build_submission as b
 manifest,genes,targets,_,base=load();effects=[];changes=[];inputs=[];routing={};pdir=A/'production_inputs';pdir.mkdir(exist_ok=True);assert not (A/'reference_model.json').exists()
 for ci,c in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(ROOT,c,genes,manifest);depth,mean,bulk,_=moments(x);q=mean*1e6;pred,traces=m.predict(q,base[ci],bulk,());e=pred[p['candidate']];effects.append(e);routing[c]=traces[p['candidate']]
  path=pdir/(c+'_query.npz');np.savez_compressed(path,query=q,query_mean=q,baseline=base[ci],genes=np.asarray(genes,str),bulk=bulk,library=library)
  inputs.append(dict(context=c,file=f'context_{c}.h5ad',sha256=sha(ROOT/f'context_{c}.h5ad'),query_sha256=sha(path),control_cells=x.shape[0]))
  for ti,t in enumerate(targets):changes.append(dict(context=c,target_gene=str(t),sum_absolute_LFC_correction=float(abs(e[ti]-base[ci,ti]).sum()),changed_genes=int(np.count_nonzero(e[ti]!=base[ci,ti]))))
  del x;gc.collect();print('V66 production',c,p['candidate'],flush=True)
 path=A/'reference_model_arrays.npz';np.savez_compressed(path,genes=np.asarray(genes,str),targets=np.asarray(targets,str),effects=np.asarray(effects,np.float32),parent_effects=base)
 pd.DataFrame(changes).to_csv(A/'production_response_changes.csv',index=False)
 write(A/'reference_model.json',dict(model='V66 strength-matched original-assay response',candidate=p['candidate'],status='public_supported_candidate' if assessment['release_supported'] else 'experimental_failed_public_guards',available_public_guards_passed=assessment['release_supported'],confirmed_improvement=False,emitter='residual10',model_arrays_sha256=sha(path),scientific_plan_sha256=sha(A/'plan.json'),production_plan_sha256=sha(A/'production_plan.json'),count_assessment_sha256=sha(A/'count_assessment.json'),context_original_inputs=inputs,dependencies=p['dependencies'],parent_V24_arrays_sha256=sha(ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),guard_traces=routing,D3_opened=False,uploaded=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit']);p.add_argument('--experimental',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  else:fit(a.experimental)
