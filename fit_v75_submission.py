"""Fit V75 production effects only after mean selection and paired count assessment."""
import gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v75_assay_shrink as m
from v57_pie import read,write,sha

def main():
 m.check();selected=read(m.A/'selection.json')['selected'];assert selected
 assessment=read(m.A/'count_assessment.json')
 # Even an experimental release must improve generated profiles and beat shuffle.
 assert all(assessment['checks'][k] for k in ('profile','worst_study','pds','beats_shuffle')),'No profile/PDS-supported production release'
 repeat=read(m.A/'count_assessment_repeat.json')
 assert all(repeat['checks'][k] for k in ('profile','worst_study','pds','beats_shuffle')),'Repeat critical gates fail'
 assert read(m.A/'count_audit.json')['passed'] and read(m.A/'count_audit_repeat.json')['passed']
 assert read(m.A/'validation_complete.json')['selection_unchanged']
 assert read(m.A/'boundary_audit.json')['passed']
 assert read(m.A/'broader_assessment.json')['passes'],'Broader response diagnostic blocks production'
 assert not (m.A/'reference_model.json').exists()
 files=['v75_assay_shrink.py','v70_count_teacher.py','v62_conserved_response.py','fit_v7_response.py','fit_v75_submission.py','build_v75_submission.py','package_v75_submission.py','validate_v75_submission.py','verify_v75_effects.py','v16_methods.py','v63_identity_guard.py']
 dependencies={n:sha(m.ROOT/n) for n in files}
 parent=m.ROOT/'experiments/v66_assay_strength/analysis';record=read(parent/'reference_model.json')
 assert sha(parent/'reference_model_arrays.npz')==record['model_arrays_sha256']
 with np.load(parent/'reference_model_arrays.npz',allow_pickle=False) as z:base=z['effects'].copy();genes=z['genes'].copy();targets=z['targets'].copy()
 inputs=[];effects=[];changes=[];pdir=m.A/'production_inputs';pdir.mkdir(exist_ok=True)
 plan=dict(candidate=selected,dependencies=dependencies,parent_model_sha256=sha(parent/'reference_model.json'),selection_sha256=sha(m.A/'selection.json'),count_assessment_sha256=sha(m.A/'count_assessment.json'),automatic_upload=False,D3_opened=False)
 write(m.A/'production_plan.json',plan)
 for ci,c in enumerate(('A','B','C')):
  old=parent/'production_inputs'/(c+'_query.npz');original=record['context_original_inputs'][ci];assert sha(old)==original['query_sha256']
  with np.load(old,allow_pickle=False) as z:q=z['query'].copy();bulk=z['bulk'].copy();library=z['library'].copy()
  e=m.predict(base[ci],bulk,())[selected];effects.append(e)
  path=pdir/(c+'_query.npz');np.savez_compressed(path,query=q,baseline=base[ci],bulk=bulk,library=library,genes=genes)
  inputs.append(dict(original,query_sha256=sha(path)))
  for ti,t in enumerate(targets):changes.append(dict(context=c,target_gene=str(t),sum_absolute_LFC_correction=float(abs(e[ti]-base[ci,ti]).sum()),changed_genes=int(np.count_nonzero(e[ti]!=base[ci,ti]))))
  print('V75 fitted production',c,flush=True);gc.collect()
 import pandas as pd
 pd.DataFrame(changes).to_csv(m.A/'production_response_changes.csv',index=False)
 path=m.A/'reference_model_arrays.npz';np.savez_compressed(path,genes=genes,targets=targets,effects=np.asarray(effects,np.float32),parent_effects=base)
 write(m.A/'reference_model.json',dict(model='V75 raw-count response teacher',candidate=selected,status='public_supported_candidate' if (assessment['release_supported'] and repeat['release_supported']) else 'experimental_incomplete_public_guards',available_public_guards_passed=(assessment['release_supported'] and repeat['release_supported']),confirmed_improvement=False,emitter='residual10',model_arrays_sha256=sha(path),scientific_plan_sha256=sha(m.A/'plan.json'),production_plan_sha256=sha(m.A/'production_plan.json'),count_assessment_sha256=sha(m.A/'count_assessment.json'),context_original_inputs=inputs,dependencies=dependencies,parent_V66_arrays_sha256=sha(parent/'reference_model_arrays.npz'),D3_opened=False,uploaded=False))

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
