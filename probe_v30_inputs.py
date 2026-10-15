"""Metadata-only audit for the new full-target mean-response pilot."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
D=ROOT/'experiments/v30_structured_response/analysis';D.mkdir(parents=True,exist_ok=True)
records=[]
for rel in ['reference/K562_gwps_raw_bulk_01.h5ad','reference/v13_resources/rpe1_raw_bulk_01.h5ad']:
 a=ad.read_h5ad(ROOT/rel,backed='r'); item={'path':rel,'shape':list(a.shape),'obs_columns':list(a.obs),'var_columns':list(a.var),'uns_keys':list(a.uns),'layers':list(a.layers),'obs_preview':a.obs.head(3).astype(str).to_dict('records'),'var_preview':a.var.head(3).astype(str).to_dict('records')}
 for col in a.obs:
  if any(t in col.lower() for t in ['gene','target','pert','guide','count','cell','expression']):
   v=a.obs[col];item[col]={'unique':int(v.nunique()),'examples':v.astype(str).unique()[:8].tolist()}
 x=a.X[:5,:]; x=x.toarray() if hasattr(x,'toarray') else np.asarray(x);item['first5_stats']={'min':float(x.min()),'max':float(x.max()),'sums':x.sum(1).tolist(),'integer_fraction':float(np.isclose(x,np.rint(x)).mean())};a.file.close();records.append(item)
p=ROOT/'reference/v18_pie/sources/string_space/meta.json';m=json.loads(p.read_text());records.append({'string_meta_keys':list(m),'string_meta_preview':{k:(v[:4] if isinstance(v,list) else v) for k,v in m.items()},'embedding_shape':list(np.load(p.with_name('embeddings.npy'),mmap_mode='r').shape)})
for rel in ['reference/v7_sources/HCT116/selected_cells.csv','reference/v7_sources/HEK293T/selected_cells.csv','reference/v5_h1/selected_cells.csv','reference/v22_cd4/selected_cells.csv','reference/v22_cd4_d2/selected_cells.csv']:
 p=ROOT/rel
 if p.exists():
  f=pd.read_csv(p,nrows=3);records.append({'metadata':rel,'columns':list(f),'first':f.astype(str).to_dict('records')})
records=json.loads(json.dumps(records,default=str),parse_constant=lambda value:None)
(D/'input_probe.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n');print(json.dumps(records[:2],indent=2));print('STRING',records[2]['embedding_shape'],records[2]['string_meta_keys'])
