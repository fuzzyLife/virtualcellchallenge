"""Metadata-only feasibility check for a fresh CD4 stimulation holdout."""
import json
from pathlib import Path
import h5py,pandas as pd
pd.options.future.infer_string=False
from probe_v22_cd4_power import Remote
from prepare_reference_controls import column
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v24_cd4_response/analysis'
def main():
 A.mkdir(parents=True,exist_ok=True);record=json.loads((ROOT/'experiments/v22_state_response/analysis/cipher_zenodo_record.json').read_text());targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy(str)
 for name in ['Marson2025_D2_Stim48hr_filtered.h5ad','Marson2025_D2_Stim8hr_filtered.h5ad']:
  item=next(v for v in record['files'] if v['key']==name);rr=Remote(item)
  with h5py.File(rr,'r') as h:
   values=column(h['obs']['perturbation']);counts=pd.Series(values).value_counts();over={str(t):int(counts.get(t,0)) for t in targets};out={'file':item,'competition_targets_present':sum(v>0 for v in over.values()),'targets400':sum(v>=400 for v in over.values()),'targets200':sum(v>=200 for v in over.values()),'counts':over,'curated_source_perturbation_categories':len(counts),'metadata_range_bytes':rr.used,'treatment_expression_scored':False}
  b.write_json(A/(name+'_power.json'),out);print(name,'present',out['competition_targets_present'],'400-cell targets',out['targets400'],flush=True)
if __name__=='__main__':main()
