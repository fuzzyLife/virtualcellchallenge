"""Inspect source label namespaces before interpreting zero literal target overlap."""
from pathlib import Path
import json,numpy as np,h5py
from prepare_reference_controls import column
from probe_v25_sources import Remote,ROOT,A

def main():
 record=json.loads((ROOT/'experiments/v24_cd4_response/analysis/cipher_zenodo_record.json').read_text());out=[]
 for item in record['files']:
  if not item['key'].startswith('kaden25'):continue
  remote=Remote(item)
  with h5py.File(remote,'r') as h:
   labels=np.unique(column(h['obs']['perturbation'])).astype(str);var=h['var'];index=var.attrs.get('_index','_index');symbols=column(var[index]).astype(str);row={'file':item['key'],'label_examples':labels[:80].tolist(),'var_columns':list(var),'gene_examples':symbols[:20].tolist(),'control_labels':[v for v in labels if 'control' in v.lower() or 'non-target' in v.lower()],'perturbation_type':'CRISPRa per primary Kaden/Southard paper; not direct CRISPRi labels'};out.append(row);print(json.dumps(row),flush=True)
 (A/'label_namespace_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
if __name__=='__main__':main()
