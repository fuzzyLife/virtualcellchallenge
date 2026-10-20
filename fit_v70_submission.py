"""Fit V70 production effects only after mean selection and paired count assessment."""
import gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as m
from v57_pie import read,write,sha

def main():
 m.check();selected=read(m.A/'selection.json')['selected'];assert selected
 assessment=read(m.A/'count_assessment.json')
 # Even an experimental release must improve generated profiles and beat shuffle.
 assert assessment['checks']['profile'] and assessment['checks']['worst_study'] and assessment['checks']['beats_shuffle'],'No profile-supported production release'
 assert read(m.A/'boundary_audit.json')['passed']
 assert read(m.A/'broader_assessment.json')['passes'],'Broader response diagnostic blocks production'
 assert not (m.A/'reference_model.json').exists()
 files=['v70_count_teacher.py','fit_v70_submission.py','build_v70_submission.py','package_v70_submission.py','validate_v70_submission.py','verify_v70_effects.py','v16_methods.py','v63_identity_guard.py']
 dependencies={n:sha(m.ROOT/n) for n in files}
 for source in m.SOURCES:
  n=f'experiments/v70_count_teacher/analysis/sources/{source}.npz';dependencies[n]=sha(m.ROOT/n)
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
  print('V70 fitted production',c,flush=True);gc.collect()
 import pandas as pd
 pd.DataFrame(changes).to_csv(m.A/'production_response_changes.csv',index=False)
 path=m.A/'reference_model_arrays.npz';np.savez_compressed(path,genes=genes,targets=targets,effects=np.asarray(effects,np.float32),parent_effects=base)
 write(m.A/'reference_model.json',dict(model='V70 raw-count response teacher',candidate=selected,status='public_supported_candidate' if assessment['release_supported'] else 'experimental_incomplete_public_guards',available_public_guards_passed=assessment['release_supported'],confirmed_improvement=False,emitter='residual10',model_arrays_sha256=sha(path),scientific_plan_sha256=sha(m.A/'plan.json'),production_plan_sha256=sha(m.A/'production_plan.json'),count_assessment_sha256=sha(m.A/'count_assessment.json'),context_original_inputs=inputs,dependencies=dependencies,parent_V66_arrays_sha256=sha(parent/'reference_model_arrays.npz'),D3_opened=False,uploaded=False))

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
