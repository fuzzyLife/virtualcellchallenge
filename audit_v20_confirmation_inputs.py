"""Audit frozen-choice reserved prompt provenance after Stack inference, without fitting."""
from pathlib import Path
import json,csv,hashlib
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def main():
 plan=json.loads((A/'plan.json').read_text());prep=json.loads((A/'confirmation_input_preparation.json').read_text());assert prep['complete'] and prep['frozen_source_choices_reused'] and prep['original_raw_counts_and_same_batch_pairing_reconstructed'] and not prep['query_treatments_used'];assert prep['plan_sha256']==sha(A/'plan.json') and prep['script_sha256']==sha(ROOT/'prepare_v20_confirmation_inputs.py')
 for item in prep['input_files']:assert sha(ROOT/item['path'])==item['sha256']
 with (ROOT/'pert_counts.csv').open() as stream:targets=[r[0] for r in list(csv.reader(stream))[1:]]
 checks=[]
 for context in plan['contexts']:
  family=context['family'];rec=json.loads((A/'inputs'/(family+'.json')).read_text());p=A/'experts'/(family+'_confirmation.json');expert=json.loads(p.read_text());assert sha(p.with_suffix('.npz'))==expert['sha256'] and expert['plan_sha256']==sha(A/'plan.json');assert expert['query_input_sha256']==rec['query_arrays_sha256'];expected={c['target'] for c in context['cells'] if c['split']=='confirmation'};assert {r['target'] for r in expert['records']}==expected
  for row in expert['records']:
   ti=targets.index(row['target']);source=rec['source_choices'][ti];assert row['source_family']==source and source!=family and row['input_count']==rec['matched_prompt_sizes'][ti]
   if not source:continue
   assert {v['kind'] for v in row['pairs']}=={'real','shuffled'}
   for pair in row['pairs']:
    target=targets[ti if pair['kind']=='real' else rec['shuffled_target_indices'][ti]];side=A/'inputs/sources'/(source+'_'+target+'.json');info=json.loads(side.read_text());assert pair['source_target']==target and pair['source_input_sha256']==sha(side);assert pair['seed']==20261020+plan['family_order'].index(family)*10000+ti;assert info['batches_match_by_row'] and info['records']['treated']['unique_original_cells']==info['source_input_rows'] and row['input_count']<=info['source_input_rows'];assert all(len(pair[k])==64 for k in ['control_mean_array_sha256','treated_mean_array_sha256'])
   checks.append({'family':family,'target':row['target'],'source_family':source,'matched_prompt_size':row['input_count']})
 output={'passed':True,'all_reserved_prompt_files_verified':len(prep['input_files']),'paired_target_groups':len(checks),'held_family_exclusion_and_exact_frozen_choices':True,'query_treatment_inputs_used':False,'checks':checks,'plan_sha256':sha(A/'plan.json'),'preparation_receipt_sha256':sha(A/'confirmation_input_preparation.json'),'auditor_sha256':sha(Path(__file__))}
 (A/'confirmation_source_input_independent_audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8');print('Reserved source input audit passed',len(checks),'groups',flush=True)
if __name__=='__main__':main()
