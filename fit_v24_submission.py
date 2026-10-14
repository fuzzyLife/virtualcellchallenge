"""Freeze a confirmed V24 response expert on exact published V16 ABC effects."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
import build_submission as b
from v24_methods import ROOT,OUT,A,load_plan,make_teacher,predict,write
from fit_v7_response import read_cache
from v6_methods import moments

def main():
 from v24_runtime import require_science
 require_science()
 plan=load_plan();selection=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());audit=json.loads((A/'independent_experiment_audit.json').read_text());assert selection['candidate']=='cdf50_de' and confirmation['candidate']==selection['candidate'] and confirmation['new_candidate_passes'] and audit['passed'];assert audit['full_axis_replayed_cases']>=4
 assert not (A/'reference_model.json').exists();data=read_cache();manifest,genes,targets=b.read_spec(ROOT);np.testing.assert_array_equal(genes,data['genes']);np.testing.assert_array_equal(targets,data['targets']);teacher=make_teacher(True);effects=[];inputs=[];changes=[];routing={};folder=A/'production_inputs';folder.mkdir(exist_ok=True)
 parent=ROOT/'experiments/v16_control_preserving/analysis';parent_model=json.loads((parent/'reference_model.json').read_text());assert b.sha256(parent/'reference_model_arrays.npz')==parent_model['model_arrays_sha256']
 with np.load(parent/'reference_model_arrays.npz',allow_pickle=False) as z:baseline=z['effects'].copy();np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
 for ci,context in enumerate(manifest['contexts']):
  x,*_=b.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6;np.savez_compressed(folder/(context+'_query.npz'),genes=genes,query_mean=query);pred,diag=predict(data,teacher,query,production=True,baseline=baseline[ci]);value=pred[selection['candidate']];assert abs(value-baseline[ci]).max()<=.250001;effects.append(value);routing[context]=diag;inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':b.sha256(ROOT/f'context_{context}.h5ad'),'query_sha256':b.sha256(folder/(context+'_query.npz'))})
  for ti,t in enumerate(targets):changes.append({'context':context,'target_gene':str(t),'sum_absolute_LFC_correction':float(abs(value[ti]-baseline[ci,ti]).sum()),'max_absolute_LFC_correction':float(abs(value[ti]-baseline[ci,ti]).max()),'changed_genes':int(np.count_nonzero(value[ti]!=baseline[ci,ti]))})
  del x;gc.collect();print('V24 production response fitted',context,diag,flush=True)
 np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects,dtype=np.float32));pd.DataFrame(changes).to_csv(A/'production_response_changes.csv',index=False)
 deps=set(plan['dependencies'])|{'v24_runtime.py','fit_v24_submission.py','build_v24_submission.py','package_v24_submission.py','validate_v24_submission.py','unpack_v24_official_payload.py','audit_v24_experiment.py','audit_v24_cd4_inputs.py'}
 write(A/'reference_model.json',{'model':'V24 context-matched donor-supported CD4 mean response','candidate':selection['candidate'],'emitter':'residual10','routing':routing,'confirmed_improvement':False,'available_public_guards_passed':True,'status':'exploratory_public_condition_transfer_verification','runtime_sha256':b.sha256(A/'runtime.json'),'model_arrays_sha256':b.sha256(A/'reference_model_arrays.npz'),'context_original_inputs':inputs,'teacher_all_sha256':b.sha256(A/'models/teacher_all.npz'),'teacher_all_record_sha256':b.sha256(A/'models/teacher_all.json'),'parent_V16_arrays_sha256':b.sha256(parent/'reference_model_arrays.npz'),'plan_sha256':b.sha256(A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'independent_experiment_audit_sha256':b.sha256(A/'independent_experiment_audit.json'),'dependencies':{n:b.sha256(ROOT/n) for n in sorted(deps)},'scope':'Development calibration overlaps public Rest summaries; condition reserve excludes Stim48hr response rows, but upstream nuisance preprocessing may overlap. All three public conditions used only in standard final refit. Private controls supply compatibility, never hidden treatment labels. Missing development NMAE prevents complete six-metric verification; no hidden leaderboard guarantee.'})
 print('V24 confirmed production response frozen',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
