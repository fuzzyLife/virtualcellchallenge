"""Supplemental checks of corrected flow implementation and count contracts."""
import json
from pathlib import Path
import v32_methods as m
D=m.ROOT/'experiments/v34_residual_flow/analysis'
plan=json.loads((D/'plan.json').read_text())
for n,h in plan['source_inputs'].items():assert m.sha(m.ROOT/n)==h,n
amend=json.loads((D/'integration_correction_amendment.json').read_text());assert m.sha(m.ROOT/'v34_corrected_methods.py')==amend['corrected_source_sha256'];assert m.sha(m.ROOT/'run_v34_corrected.py')==amend['wrapper_sha256'];rec=json.loads((D/'count_projection_recovery.json').read_text());assert m.sha(m.ROOT/'run_v34_projected_counts.py')==rec['recovery_code_sha256']
rows=[]
for stage in ['null','treatment']:
 data=json.loads((D/f'{stage}_results.json').read_text())
 for r in data['records']:
  if r['arm']=='flow_shape25':
   assert r['gene_totals_match_mean_reference'] and r['mean_reference_L1']<=.0005;rows.append({'stage':stage,'family':r['family'],'target':r['target'],'CPM_L1':r['mean_reference_L1']})
assert len(rows)==25
m.write(D/'count_contract_audit.json',{'passed':True,'shape_populations':len(rows),'maximum_CPM_L1':max(v['CPM_L1'] for v in rows),'records':rows,'checks':'Frozen source-cache and teacher hashes; corrected integration/count code hashes; all25 reported exact-total and strict CPM contracts checked. Actual integer libraries/totals were asserted during every generated population. Counts are reproducible from seeds/models but not independently regenerated in this supplement.','auditor_sha256':m.sha(m.ROOT/'audit_v34_contracts.py')});print('V34 count contracts checked',len(rows),flush=True)