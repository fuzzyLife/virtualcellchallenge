"""Fetch only newly required CD4 DE rows using the original pinned source identities."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json, numpy as np,pandas as pd,h5py
import extract_cd4_reference as source
import build_submission as b
R=Path(__file__).resolve().parent;O=R/'reference/final_v91/cd4';A=R/'experiments/final_v91/analysis'
def main():
 O.mkdir(parents=True,exist_ok=True);A.mkdir(parents=True,exist_ok=True)
 new=np.array(sorted(set(pd.read_csv(R/'final/pert_counts.csv').target_gene)-set(pd.read_csv(R/'pert_counts.csv').target_gene)))
 records=[]
 for kind,identity in source.OBJECTS.items():
  reader=source.MetadataReader(identity)
  with h5py.File(reader,'r') as f:
   groups={'merged':f} if kind=='merged' else {n:f['mod'][n] for n in ['CE0006864_CE0008162','CE0008678_CE0010866']}
   for name,group in groups.items():
    dest=O/(name+'_new.npz');receipt=dest.with_suffix('.json')
    if not receipt.exists():
     rec=source.extract(group,identity,new,dest);b.write_json(receipt,rec)
    rec=json.loads(receipt.read_text());assert b.sha256(dest)==rec['sha256'];records.append(rec)
    oldmeta=pd.read_csv(R/'reference/cd4_probe'/(name+'.csv'));newmeta=pd.read_csv(dest.with_suffix('.csv'))
    with np.load(R/'reference/cd4_probe'/(name+'.npz'),allow_pickle=False) as z:old={k:z[k].copy() for k in z.files}
    with np.load(dest,allow_pickle=False) as z:extra={k:z[k].copy() for k in z.files}
    np.testing.assert_array_equal(old['genes'],extra['genes']);meta=pd.concat([oldmeta,newmeta],ignore_index=True)
    assert meta.source_row.is_unique;order=np.argsort(meta.source_row.to_numpy());meta=meta.iloc[order].reset_index(drop=True);arrays={k:old[k] if k=='genes' else np.concatenate([old[k],extra[k]],axis=0)[order] for k in old}
    np.testing.assert_array_equal(meta.source_row.to_numpy(),arrays['source_rows']);np.savez_compressed(O/(name+'.npz'),**arrays);meta.to_csv(O/(name+'.csv'),index=False)
    print('Final CD4 source assembled',name,len(meta),flush=True)
  reader.close()
 b.write_json(A/'cd4_collection.json',dict(new_targets=new.tolist(),new_receipts=records,combined={n:b.sha256(O/n) for n in ['merged.npz','merged.csv','CE0006864_CE0008162.npz','CE0006864_CE0008162.csv','CE0008678_CE0010866.npz','CE0008678_CE0010866.csv']},collector_sha256=b.sha256(Path(__file__))))
if __name__=='__main__':main()
