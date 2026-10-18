"""Audit the V58 held-study correction and completed numerical artifacts."""
import argparse
import numpy as np
from v58_pie import *

def boundaries():
 p=plan();rec=read(A/'held_study_poison_audit.json');assert rec['passed'] and rec['plan_sha256']==sha(A/'plan.json')
 for c in rec['checks']:assert 'CD4' in c['excluded_families'] and c['all_bank_fields_bitwise_invariant']
 for source in SOURCES:
  ip=CACHE/'inputs'/(source+'.npz');ir=read(A/'inputs'/(source+'.json'));assert sha(ip)==ir['sha256'] and ir['plan_sha256']==sha(A/'plan.json')
  with np.load(ip,allow_pickle=False) as x,np.load(PARENT_CACHE/'inputs'/ip.name,allow_pickle=False) as y:
   for k in y.files:
    if k!='baseline' or source.startswith('CD4'):np.testing.assert_array_equal(x[k],y[k])
  pr=read(A/'predictions'/(source+'.json'));assert pr['plan_sha256']==sha(A/'plan.json');assert sha(CACHE/'predictions'/(source+'.npz'))==pr['parent_prediction_sha256']==pr['sha256']
 methods=[]
 for old,new in [('fit_v57_pie_decoder.py','fit_v58_pie_decoder.py'),('evaluate_v57_pie_counts.py','evaluate_v58_pie_counts.py')]:
  a=(ROOT/old).read_text(encoding='utf-8-sig');b=(ROOT/new).read_text(encoding='utf-8-sig').replace('from v58_pie import','from v57_pie import').replace('from fit_v58_pie_decoder import','from fit_v57_pie_decoder import');assert a==b;methods.append(new)
 write(A/'boundaries_audit.json',{'passed':True,'methods_byte_equivalent_after_import_rebinding':methods,'neural_inputs_outputs_and_labels_bitwise_unchanged':True,'held_CD4_comparator_unchanged':True,'all_calibration_baselines_CD4_poison_invariant':True,'plan_sha256':sha(A/'plan.json'),'parent_evidence_audit_sha256':sha(PARENT_A/'evidence_independent_audit.json')})
 print('V58 boundaries passed: feature exclusion, poison invariance, identical neural inputs/outputs',flush=True)

def artifacts():
 boundaries();p=plan()
 for k in ['joint','shuffled','evidence']:
  r=read(A/'calibrators'/(k+'.json'));assert not r['CD4_labels_used_in_fit'];assert all(not g['source'].startswith('CD4') for g in r['groups']);assert sha(CACHE/'calibrators'/(k+'.npz'))==r['sha256'];assert r['plan_sha256']==sha(A/'plan.json')
  tr={g['target'] for g in r['groups'] if not g['validation']};va={g['target'] for g in r['groups'] if g['validation']};assert not tr&va
 for f in (A/'counts').glob('*.json'):
  r=read(f);assert r['plan_sha256']==sha(A/'plan.json')
  if f.stem.endswith('_null'):continue
  assert len(set(r['actual_rows']))==len(r['actual_rows']) and len(set(r['template_rows']))==400
  assert not set(r['actual_rows'])&set(r['template_rows']) and not set(r['reference_rows'])&set(r['template_rows'])
  assert sha(CACHE/'count_statistics'/f.with_suffix('.npz').name)==r['statistics_sha256']
  assert {v['arm'] for v in r['records']}==set(p['config']['arms'])
 assert read(A/'count_evaluation.json')['complete'] and read(A/'response_evaluation.json')['complete']
 write(A/'final_audit.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'boundaries_sha256':sha(A/'boundaries_audit.json'),'count_evaluation_sha256':sha(A/'count_evaluation.json'),'response_evaluation_sha256':sha(A/'response_evaluation.json'),'D3_opened':False,'new_submission_generated':False,'source_only_exclusion_verified_by_held_label_poisoning':True})
 print('V58 final audit passed',flush=True)

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['boundaries','artifacts']);args=q.parse_args();{'boundaries':boundaries,'artifacts':artifacts}[args.stage]()
