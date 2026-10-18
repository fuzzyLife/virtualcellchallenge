"""Seal the V58 exclusion repair, poison-test it and reuse unchanged PIE outputs."""
import argparse,copy,shutil,gc
import numpy as np
from v58_pie import *

def freeze():
 assert not (A/'plan.json').exists()
 p=copy.deepcopy(read(PARENT_A/'plan.json'));p['version']=58;p['parent_plan_sha256']=sha(PARENT_A/'plan.json')
 p['correction']='Exclude held CD4 responses from calibration baseline features as well as from fitted labels. Source targets still exclude their own recipient study. No model architecture, alpha grid, target split, cohort, seed, feature formula or decoder change.'
 for n in ['v58_pie.py','prepare_v58_pie.py','fit_v58_pie_decoder.py','evaluate_v58_pie_counts.py','audit_v58_pie.py','run_v58_cached_stage.py']:
  p['dependencies'][n]=sha(ROOT/n)
 p['parent_notice_sha256']=sha(PARENT_A/'holdout_correction_notice.json')
 p['calibration_feature_exclusion']='H1: H1+CD4; HCT116/HEK293T: both Orion lines+CD4; CD4: all CD4. No held-study response bank values may alter any prepared baseline.'
 write(A/'plan.json',p);print('V58 frozen',sha(A/'plan.json'),flush=True)

def prepare():
 import v38_response_hurdle as hist
 from v16_methods import source_arrays
 from v8_methods import gated
 p=plan();d=hist.data();hold=d['row_families']=='CD4';poison=dict(d)
 for k,value in [('truth',17.),('errors',19.)]:
  poison[k]=d[k].copy();poison[k][hold]=value
 banks={};checks=[]
 for group in ['H1','Orion','CD4']:
  excluded=hist.legacy(tuple(sorted({group,'CD4'})));a=source_arrays(d,excluded);b=source_arrays(poison,excluded)
  for k in a:np.testing.assert_array_equal(a[k],b[k])
  banks[group]=a;checks.append({'recipient_group':group,'excluded_families':excluded,'poisoned_CD4_response_rows':int(hold.sum()),'all_bank_fields_bitwise_invariant':True})
  print('Held-study poison guard passed',group,flush=True)
 del poison;gc.collect()
 for source in SOURCES:
  original=PARENT_CACHE/'inputs'/(source+'.npz');receipt=read(PARENT_A/'inputs'/(source+'.json'));assert sha(original)==receipt['sha256']
  with np.load(original,allow_pickle=False) as z:values={k:z[k].copy() for k in z.files}
  group='CD4' if source.startswith('CD4') else ('Orion' if source in ['HCT116','HEK293T'] else 'H1');bank=banks[group]
  base=np.asarray([gated(v.copy(),np.nan_to_num(values['query_cpm']),int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32)
  if source.startswith('CD4'):np.testing.assert_array_equal(base,values['baseline'])
  changed=int(np.count_nonzero(base!=values['baseline']));values['baseline']=base
  dest=CACHE/'inputs'/(source+'.npz');dest.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(dest,**values)
  # The frozen neural input excludes baseline, labels and validation statistics.
  with np.load(original,allow_pickle=False) as z:
   for k in ['genes','targets','ctrl_mean','support','measured','query_cpm','mean','rank','valid','delta']:np.testing.assert_array_equal(values[k],z[k])
  receipt.update({'sha256':sha(dest),'plan_sha256':sha(A/'plan.json'),'parent_input_sha256':sha(original),'parent_plan_sha256':sha(PARENT_A/'plan.json'),'baseline_changed_coordinates':changed,'held_CD4_excluded_from_baseline_features':True,'original_PIE_input_fields_identical':True})
  write(A/'inputs'/(source+'.json'),receipt)
  pp=PARENT_CACHE/'predictions'/(source+'.npz');rr=read(PARENT_A/'predictions'/(source+'.json'));assert sha(pp)==rr['sha256']
  out=CACHE/'predictions'/pp.name;out.parent.mkdir(exist_ok=True);shutil.copyfile(pp,out);assert sha(out)==rr['sha256']
  rr.update({'plan_sha256':sha(A/'plan.json'),'input_sha256':sha(dest),'parent_prediction_sha256':sha(pp),'PIE_forwards_reused_unchanged':True,'reuse_reason':'Only the baseline field changed; it is not read by PIE inference. All control-only neural inputs and labels were checked bitwise equal.'})
  write(A/'predictions'/(source+'.json'),rr)
  shutil.copyfile(PARENT_A/'predictions'/(source+'_evidence_ablation.json'),A/'predictions'/(source+'_evidence_ablation.json'))
  print('V58 prepared',source,'baseline coordinates changed',changed,flush=True)
 write(A/'held_study_poison_audit.json',{'passed':True,'checks':checks,'plan_sha256':sha(A/'plan.json'),'all_CD4_comparator_baselines_unchanged':True,'all_PIE_inputs_and_outputs_reused_bitwise':True,'labels_or_control_rows_changed':False,'D3_opened':False})

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['freeze','prepare']);z=q.parse_args();{'freeze':freeze,'prepare':prepare}[z.stage]()
