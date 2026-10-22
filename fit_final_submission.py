"""Freeze and fit a final-panel adaptation of the best published validation method."""
from pathlib import Path
import argparse,gc,sys,importlib.metadata
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
import final_response_methods as f
from v57_pie import read,write,sha
R=f.R;A=f.A
FILES=['run_final_pipeline.py','test_final_cli.py','audit_final_model.py','inspect_final_inputs.py','collect_final_cd4.py','collect_final_orion.py','collect_final_h1.py','build_final_sources.py','audit_final_replay.py','test_final_count_replay.py','diagnose_final_contexts.py','audit_final_sources.py','final_response_methods.py','fit_final_submission.py','build_final_submission.py','package_final_submission.py','validate_final_submission.py','v24_methods.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py','v8_methods.py','v6_methods.py','v62_conserved_response.py','v63_identity_guard.py','v86_joint_posterior.py','v88_mixture_quantiles.py','build_submission.py','extract_cd4_reference.py','prepare_v7_sources.py','audit_v7_source_offsets.py']
def freeze():
 assert not (A/'plan.json').exists();inp=read(A/'input_audit.json');assert inp['contexts']==['D','E','F']
 for c in ['A','B','C']:assert read(A/f'validation_replay_{c}.json')['passed']
 assert read(A/'count_runtime_replay.json')['passed'];assert read(A/'source_audit.json')['passed'];assert read(A/'cli_final_smoke.json')['passed']
 receipts=['input_audit.json','source_provenance.json','cd4_collection.json','H1_new_source.json','HCT116_new_source.json','HEK293T_new_source.json','validation_replay_A.json','validation_replay_B.json','validation_replay_C.json','count_runtime_replay.json','source_audit.json','cli_final_smoke.json']
 write(A/'plan.json',dict(version='final_v91',partition='test',panel='vcc2026-test-1',contexts=inp['contexts'],method='V24 cdf50_de → V66 assay75 with 80% identity guard → V90/V86 interval25 with the same 80% guard → unchanged residual10 dual-moment count emitter.',selection='Use best published validation method V90 (0.15214377976638518). No final treatment labels, final scores or strength search. Sources are extended to the new final panel; source priors and CD4 teacher are refit on genuine public evidence. Missing observations remain masked.',changes='New D/E/F controls, new final target list, rebuilt source-target rows and corresponding source-only priors/teacher, newly fitted source-panel centering and rank axes. Numerical runtime is recorded and old validation response replay required within absolute 1e-5 LFC. No claim of bitwise equality across runtimes or of comparable final scores.',template_seed=20261011,own_target_factor=.1,emitter='residual10',finite_prior_EM_updates=100,mixture_fit_coordinates=20000,interval_normal_coverage=.8,interval_projection_strength=.25,maximum_interval_LFC_correction=.125,automatic_upload=False,score_visibility='blind',independent_final_accuracy_validation=False,code_sha256={n:sha(R/n) for n in FILES},source_receipt_sha256={str((A/n).relative_to(R)):sha(A/n) for n in receipts},published_validation_receipt_sha256=sha(R/'experiments/v90_interval_diagnostic/analysis/published_result.json'),runtime=dict(python=sys.version,packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','pandas','anndata','h5py','threadpoolctl','vcc-cli']})))
 print('Final V91 method frozen before final inference',flush=True)
def prior(shuffled=False):
 f.check();f.posterior().fit('production',shuffled)
def fit():
 f.check();assert not (A/'reference_model.json').exists();m,g,t=b.read_spec(R/'final');effects=[];parents=[];trace={};inputs=[]
 for c in m['contexts']:
  with np.load(A/(c+'_controls.npz'),allow_pickle=False) as z:q=z['query'].copy();bulk=z['bulk'].copy();lib=z['library'].copy()
  pred,tr,old=f.predict(q,bulk);effects.append(pred);parents.append(old['v66']);trace[c]=tr;np.savez_compressed(A/(c+'_response.npz'),query=q,bulk=bulk,library=lib,final=pred,**old)
  inputs.append(dict(context=c,file='final/context_'+c+'.h5ad',sha256=sha(R/'final'/('context_'+c+'.h5ad')),query_sha256=sha(A/(c+'_controls.npz')),response_sha256=sha(A/(c+'_response.npz'))));print('Final fitted response',c,flush=True);gc.collect()
 path=A/'reference_model_arrays.npz';np.savez_compressed(path,genes=g.astype(str),targets=t.astype(str),effects=np.asarray(effects,np.float32),parent_effects=np.asarray(parents,np.float32));p=read(A/'plan.json');write(A/'reference_model.json',dict(model='Final-panel transfer of published V90 method',candidate='interval25',partition='test',contexts=m['contexts'],status='blind_final_transfer_no_final_accuracy_claim',emitter='residual10',model_arrays_sha256=sha(path),plan_sha256=sha(A/'plan.json'),context_original_inputs=inputs,dependencies=p['code_sha256'],guard_traces=trace,uploaded=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prior','fit']);p.add_argument('--shuffled',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  elif a.stage=='prior':prior(a.shuffled)
  else:fit()
