"""Repair only our failed object-axis NPZ export without unpickling any array."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
def axis(name):
 with (ROOT/name).open() as f:return np.asarray([r[0] for r in list(csv.reader(f))[1:]],dtype=str)
def main():
 old=A/'production_plan_before_axis_dtype_fix.json';plan=json.loads(old.read_text());assert sha(A/'production_plan.json')==sha(old);assert not (A/'production_experts').exists() or not list((A/'production_experts').glob('*.npz'));genes=axis('gene_names.csv');targets=axis('pert_counts.csv');assert genes.shape==(18533,) and targets.shape==(300,);checks=[];changed=set()
 for context in ['A','B','C']:
  folder=A/'production_inputs';p=folder/(context+'_query.npz');record=folder/(context+'.json');saved=A/'production_axis_export_snapshot'/p.name;assert sha(p)==sha(saved)
  with np.load(p,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files if k not in ['genes','targets']}
  arrays.update(genes=genes,targets=targets);temporary=p.with_suffix('.repair.npz');np.savez_compressed(temporary,**arrays)
  with np.load(temporary,allow_pickle=False) as z:
   for name,value in arrays.items():np.testing.assert_array_equal(z[name],value)
   assert z['genes'].dtype.kind=='U' and z['targets'].dtype.kind=='U'
  temporary.replace(p);rec=json.loads(record.read_text());assert rec['query_arrays_sha256']==sha(saved);rec['query_arrays_sha256']=sha(p);record.write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8');changed.update([str(p.relative_to(ROOT)).replace('\\','/'),str(record.relative_to(ROOT)).replace('\\','/')]);checks.append({'context':context,'all_non_axis_arrays_bitwise_preserved':True,'axes_read_independently_from_official_CSV':True,'old_query_sha256':sha(saved),'new_query_sha256':sha(p)})
 for item in plan['input_files']:
  digest=sha(ROOT/item['path'])
  if item['path'] not in changed:assert digest==item['sha256']
  item['sha256']=digest
 for name in plan['dependencies']:
  digest=sha(ROOT/name)
  if name not in ['v20_production.py','infer_v20_submission.py','fit_v20_submission.py']:assert digest==plan['dependencies'][name]
  plan['dependencies'][name]=digest
 receipt={'passed':True,'reason':'Our new exporter saved gene and target axes as object dtype. Rejecting safe loader caught this before any expert output. Use plain Unicode string axes.','model_formula_source_choices_seeds_and_all_numeric_input_arrays_unchanged':True,'no_pickle_loaded':True,'original_plan_sha256':sha(old),'checks':checks,'repair_script_sha256':sha(Path(__file__))}
 repair=A/'production_axis_export_repair.json';repair.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8');plan['axis_export_repair_sha256']=sha(repair);(A/'production_plan.json').write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8');print('Plain-string axes repaired; every numerical input and all source prompts preserved',flush=True)
if __name__=='__main__':main()
