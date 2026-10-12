"""Collect pinned public metadata before downloading large expression matrices."""
import json,requests,hashlib
import numpy as np,pandas as pd
import build_submission as b
ROOT=b.ROOT; A=ROOT/'experiments/v17_learned_heterogeneity/analysis'
def main():
 A.mkdir(parents=True,exist_ok=True)
 api='https://huggingface.co/api/datasets/arcinstitute/PIE_replogle_nadig_essential'
 info=requests.get(api,timeout=60);info.raise_for_status();revision=info.json()['sha']
 url=f'https://huggingface.co/datasets/arcinstitute/PIE_replogle_nadig_essential/resolve/{revision}/preprocessed/meta.json'
 r=requests.get(url,timeout=60);r.raise_for_status();file=A/'nadig_public_metadata.json';file.write_bytes(r.content);d=r.json()
 _,genes,targets=b.read_spec(ROOT)
 print('Metadata keys',list(d),flush=True)
 vocab=list(d['pert_to_id'])
 if not vocab:
  for k,v in d.items():
   if isinstance(v,list) and k not in ['genes','contexts']:print(k,str(v)[:200],flush=True)
 overlap=sorted(set(targets)&set(vocab))
 records={'retrieved_date':'2026-10-12','dataset':'Arc PIE Replogle/Nadig essential-gene screens','revision':revision,'url':url,'sha256':b.sha256(file),'bytes':len(r.content),'metadata_keys':list(d),'perturbation_count':len(vocab),'competition_target_overlap':overlap,'gene_overlap':len(set(genes)&set(d['genes'])),'raw_counts_downloaded':False,'reason':'Inspect exact target overlap before 5.2/8.7 GB raw downloads. Summary metadata is not raw cell-level training data.'}
 b.write_json(A/'public_data_audit.json',records);print(json.dumps(records,indent=2),flush=True)
if __name__=='__main__':main()
