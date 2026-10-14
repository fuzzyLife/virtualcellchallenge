"""V21 Stack-compatible control programs; all V20 response/emitter constants retained."""
import json
import numpy as np
import build_submission as b
import v20_methods as parent
ROOT=b.ROOT;OUT=ROOT/'experiments/v21_stack_programs';A=OUT/'analysis';R=parent.R;FAMILIES=parent.FAMILIES
NAMES=['v16','stack_mean','stack_joint','shuffled_joint','reweight','amplitude']

def target_effects(*args,**kwargs):
 effects=parent.target_effects(*args,**kwargs);effects['shuffled_joint']=effects['shuffled'];return effects

def emit(template,effect,own,seed,kind,prepared,atlas,expert=None,ti=None,query_inputs=None,template_rows=None):
 if kind=='shuffled_joint':
  expert=dict(expert,program_delta=expert['shuffled_program_delta']);kind='stack_joint'
 return parent.emit(template,effect,own,seed,kind,prepared,atlas,expert,ti,query_inputs,template_rows)

def load_plan():
 plan=json.loads((A/'plan.json').read_text())
 for name,h in plan['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V21 dependency changed: '+name
 assert b.sha256(parent.A/'plan.json')==plan['parent_plan_sha256'];return plan

def load_query(family):
 file=A/'inputs'/f'{family}_query.npz';rec=json.loads((A/'inputs'/f'{family}.json').read_text());assert rec['plan_sha256']==b.sha256(A/'plan.json') and rec['query_arrays_sha256']==b.sha256(file)
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def load_expert(family,stage):
 file=A/'experts'/f'{family}_{stage}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
