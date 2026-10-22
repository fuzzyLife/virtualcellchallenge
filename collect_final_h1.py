"""Extract newly required final targets from public 2025 H1, reusing NTC cells only once."""
from pathlib import Path
import json,gc,inspect,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,h5py
from scipy import sparse
from prepare_reference_controls import RangeReader,column,index,select_stratified
from extract_v5_h1 import prefetch_v5_sparse,read_rows,SPEC
import build_submission as b
R=Path(__file__).resolve().parent;O=R/'reference/final_v91/H1';A=R/'experiments/final_v91/analysis'

# Preserve the reviewed block plan and byte checks; preload only those verified pages.
# This removes serial UNC cache reads from HDF5 access without changing any source values.
_prefetch_namespace=dict(prefetch_v5_sparse.__globals__)
_prefetch_code=inspect.getsource(prefetch_v5_sparse).replace('reader.prefetch_planned_bytes=planned','reader.prefetch_planned_bytes=planned;reader.final_planned_blocks=sorted(blocks)')
_prefetch_code=_prefetch_code.replace('missing=[b for b in sorted(blocks) if not cached_valid(b)]', 'ordered=sorted(blocks)\n with ThreadPoolExecutor(max_workers=16) as check_pool: good=list(check_pool.map(cached_valid,ordered))\n missing=[b for b,ok in zip(ordered,good) if not ok]')
exec(compile(_prefetch_code,str(R/'extract_v5_h1.py'),'exec'),_prefetch_namespace)
_prefetch=_prefetch_namespace['prefetch_v5_sparse']
def warm_prefetch(reader,x,ptr,rows):
 _prefetch(reader,x,ptr,rows)
 wanted=set(reader.final_planned_blocks)|set(reader.blocks)
 total=sum(min(reader.block_size,reader.size-k*reader.block_size) for k in wanted)
 if total>reader.budget:raise RuntimeError(f'Planned source ranges including metadata {total} exceed byte budget')
 def load(k):
  path=reader.cache/f'{k:09d}.bin';raw=path.read_bytes();expected=min(reader.block_size,reader.size-k*reader.block_size)
  assert len(raw)==expected and hashlib.sha256(raw).hexdigest()==path.with_suffix('.sha256').read_text().strip()
  return k,raw
 missing=sorted(wanted-set(reader.blocks));print('Warming verified H1 cache pages',len(missing),flush=True)
 with ThreadPoolExecutor(max_workers=16) as pool:
  for k,raw in pool.map(load,missing):reader.blocks[k]=raw
 reader.used=total;assert sum(map(len,reader.blocks.values()))==total
 print('H1 cache pages ready; same pinned bytes, bounded memory',total,flush=True)
def main():
 O.mkdir(parents=True,exist_ok=True);new=sorted(set(pd.read_csv(R/'final/pert_counts.csv').target_gene)-set(pd.read_csv(R/'pert_counts.csv').target_gene));records=[]
 for split in ['train','validation','test']:
  out=O/split;out.mkdir(exist_ok=True);receipt=out/'extraction.json'
  if receipt.exists():
   rec=json.loads(receipt.read_text());assert b.sha256(out/'counts.npz')==rec['sha256'];records.append(rec);continue
  spec=SPEC['h1-'+split];cache=R/('reference/control_profiles/ranges_H1' if split=='train' else 'reference/upstream_review/ranges_H1_'+split);reader=RangeReader(spec['url'],cache,budget=4*1024**3)
  try:
   assert reader.size==spec['bytes']
   with h5py.File(reader,'r') as f:
    meta=pd.DataFrame({k:column(f['obs'][k]) for k in ['target_gene','guide_id','batch']});meta['source_row']=np.arange(len(meta));meta['source_cell_id']=index(f['obs']);genes=index(f['var']);x=f['X'];ptr=x['indptr'][:];mask=meta.target_gene.isin(new).to_numpy()
    if not mask.any():
     b.write_json(out/'no_targets.json',dict(source=reader.identity,targets=0));continue
    for n in [16,8,4,2,1]:
     cells=select_stratified(meta,mask,['target_gene','batch'],n).sort_values('source_row').reset_index(drop=True)
     try:warm_prefetch(reader,x,ptr,cells.source_row.to_numpy());break
     except RuntimeError as e:
      if 'Planned source ranges' not in str(e):raise
      if n==1:raise
    matrix=read_rows(x,ptr,cells.source_row.to_numpy(),len(genes));cells['public_split']=split;cells.to_csv(out/'cells.csv',index=False);pd.DataFrame({'gene_name':genes}).to_csv(out/'genes.csv',index=False);sparse.save_npz(out/'counts.npz',matrix);rec=dict(split=split,source=reader.identity,cells=len(cells),targets=int(cells.target_gene.nunique()),per_target_batch=n,sha256=b.sha256(out/'counts.npz'),cells_sha256=b.sha256(out/'cells.csv'),genes_sha256=b.sha256(out/'genes.csv'),new_targets_only=True);b.write_json(receipt,rec);records.append(rec);print('Final H1 new targets',split,len(cells),rec['targets'],flush=True);del matrix;gc.collect()
  finally:reader.close()
 b.write_json(A/'H1_new_source.json',dict(splits=records,collector_sha256=b.sha256(Path(__file__)),controls='Reuse original verified H1 NTC rows once, never repeat them across splits.'))
if __name__=='__main__':main()
