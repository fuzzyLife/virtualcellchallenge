"""Verify the numerical cache amendment and production model lineage."""
import json
from pathlib import Path
import v35_methods as m

def main():
 rec=json.loads((m.P36/'normalization_cache_amendment.json').read_text());assert rec['parent_plan_sha256']==m.sha(m.A/'plan.json');assert rec['cache_code_sha256']==m.sha(m.ROOT/'run_v35_v36_cached.py');assert rec['original_methods_sha256']==m.sha(m.ROOT/'v36_methods.py');assert rec['equivalence_checks']['toy_bit_exact'] and rec['equivalence_checks']['actual_raw_subset_bit_exact'];assert rec['V36_outcomes_available'] is False
 models=[]
 for path in (m.P36/'models').glob('*.json'):
  r=json.loads(path.read_text());parent=m.A/'models'/(r['outer']+'.npz');assert m.sha(parent)==r['extra']['parent_model_sha256'];assert not set(r['extra']['population_fit_targets'])&set(r['extra']['population_validation_targets']);assert all(m.target_split(t)!=0 for t in r['extra']['population_fit_targets']);assert all(m.target_split(t)==0 for t in r['extra']['population_validation_targets']);assert r['extra']['population_manifest_sha256']==m.sha(m.P36/'populations'/(r['outer']+'.json'));models.append(r['outer'])
 assert set(m.OUTERS)<=set(models)
 m.write(m.P36/'normalization_cache_audit.json',{'passed':True,'cache_amendment_sha256':m.sha(m.P36/'normalization_cache_amendment.json'),'models_checked':models,'checks':'Bit-exact normalization checks; original frozen methods and cache hashes; all population fit/validation target boundaries and V35 parent checkpoint lineage.','auditor_sha256':m.sha(m.ROOT/'audit_v35_v36_cache.py')});print('Cache and parent-model audit passed',flush=True)
if __name__=='__main__':main()