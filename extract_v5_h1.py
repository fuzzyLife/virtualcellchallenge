"""Extract bounded public H1 validation/test responses and verify reused controls by ID/counts."""
from pathlib import Path
import argparse,json,gc,hashlib
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from prepare_reference_controls import RangeReader,column,index,select_stratified,prefetch_sparse
import build_submission as build
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'reference/v5_h1'
SPEC=json.loads((ROOT/'Virtual-Cell-Challenge-2026/sources.json').read_text())


def prefetch_v5_sparse(reader,x,indptr,rows):
 """Plan only the HDF5 storage chunks covering selected control rows."""
 from concurrent.futures import ThreadPoolExecutor
 blocks=set(reader.blocks)
 for key in ['data','indices']:
  ds=x[key];print('Storage',key,ds.shape,ds.chunks,ds.compression,flush=True)
  spans=[]
  if ds.chunks is None:
   offset=ds.id.get_offset()
   if offset is None:raise ValueError('Unknown sparse storage')
   spans=[(offset+int(indptr[r])*ds.dtype.itemsize,offset+int(indptr[r+1])*ds.dtype.itemsize) for r in rows]
  else:
   chunks=set()
   for r in rows:
    lo,hi=int(indptr[r]),int(indptr[r+1])
    chunks.update(range(lo//ds.chunks[0],(hi-1)//ds.chunks[0]+1))
   for chunk in sorted(chunks):
    info=ds.id.get_chunk_info_by_coord((chunk*ds.chunks[0],))
    if info.byte_offset is None:raise ValueError('Unallocated selected sparse chunk')
    spans.append((info.byte_offset,info.byte_offset+info.size))
  for lo,hi in spans:
   blocks.update(range(lo//reader.block_size,(hi-1)//reader.block_size+1))
 planned=sum(min(reader.block_size,reader.size-b*reader.block_size) for b in blocks)
 if planned>reader.budget:raise RuntimeError(f'Planned source ranges {planned} exceed byte budget')
 def cached_valid(block):
  path=reader.cache/f'{block:09d}.bin';sidecar=path.with_suffix('.sha256');expected=min(reader.block_size,reader.size-block*reader.block_size)
  return path.exists() and sidecar.exists() and path.stat().st_size==expected and hashlib.sha256(path.read_bytes()).hexdigest()==sidecar.read_text().strip()
 missing=[b for b in sorted(blocks) if not cached_valid(b)]
 print('Prefetch',len(missing),'blocks; bounded planned bytes',planned,flush=True)
 def fetch(block):
  lo=block*reader.block_size;hi=min(reader.size,lo+reader.block_size)
  # A separate retrying session per worker request avoids shared session mutation.
  client=requests.Session();client.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,500,502,503,504])))
  with client.get(reader.url,headers={'Range':f'bytes={lo}-{hi-1}','If-Match':reader.etag},timeout=(15,90)) as response:
   if response.status_code!=206 or response.headers.get('Content-Range')!=f'bytes {lo}-{hi-1}/{reader.size}':raise RuntimeError('Prefetch range mismatch')
   raw=response.content
  client.close()
  if len(raw)!=hi-lo:raise RuntimeError('Truncated prefetch')
  path=reader.cache/f'{block:09d}.bin';path.write_bytes(raw);path.with_suffix('.sha256').write_text(hashlib.sha256(raw).hexdigest())
  return len(raw)
 with ThreadPoolExecutor(max_workers=24) as pool:
  for i,size in enumerate(pool.map(fetch,missing)):
   reader.downloaded+=size
   if (i+1)%128==0:print('Prefetched',i+1,'/',len(missing),'blocks',flush=True)
 reader.prefetch_planned_bytes=planned


def read_rows(x,indptr,rows,n_genes):
    data=[];columns=[];pointer=[0]
    for count,row in enumerate(rows):
        lo,hi=int(indptr[row]),int(indptr[row+1]);raw=x['data'][lo:hi]
        if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any() or raw.sum()<=0:
            raise ValueError('Source does not contain nonempty raw integer counts')
        data.append(raw.astype(np.int32));columns.append(x['indices'][lo:hi].astype(np.int32));pointer.append(pointer[-1]+len(raw))
        if (count+1)%256==0:print('Extracted rows',count+1,'/',len(rows),flush=True)
    matrix=sparse.csr_matrix((np.concatenate(data),np.concatenate(columns),np.array(pointer,dtype=np.int64)),shape=(len(rows),n_genes))
    matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices();return matrix


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--budget-mib',type=int,default=4096);args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    old=ROOT/'reference/h1_response';provenance=json.loads((old/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:
        assert build.sha256(old/name)==provenance[key]
    training=sparse.load_npz(old/'selected_counts.npz');tm=pd.read_csv(old/'selected_cells.csv')
    genes=pd.read_csv(old/'gene_names.csv').gene_name.to_numpy();official=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    reader=RangeReader(SPEC['h1-train']['url'],ROOT/'reference/control_profiles/ranges_H1',budget=128*1024**2)
    with h5py.File(reader,'r') as handle:tm['source_cell_id']=index(handle['obs'])[tm.source_row.to_numpy()]
    reader.close();tm['public_split']='train'
    controls=tm[tm.target_gene=='non-targeting'];control_ids=controls.source_cell_id.to_numpy()
    control_counts=training[controls.index.to_numpy()]
    pd.DataFrame({'gene_name':genes}).to_csv(OUT/'gene_names.csv',index=False)
    records=[];all_meta=[tm];all_counts=[training]
    for split in ['validation','test']:
        directory=OUT/split;directory.mkdir(exist_ok=True);saved=directory/'extraction.json'
        if saved.exists():
            record=json.loads(saved.read_text())
            assert record['source']['url']==SPEC['h1-'+split]['url'] and record['source']['bytes']==SPEC['h1-'+split]['bytes']
            for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256')]:assert build.sha256(directory/name)==record[key]
            matrix=sparse.load_npz(directory/'selected_counts.npz');meta=pd.read_csv(directory/'selected_cells.csv')
            print('Verified cached split',split,flush=True)
        else:
            spec=SPEC['h1-'+split];cache=ROOT/'reference/upstream_review'/('ranges_H1_'+split)
            reader=RangeReader(spec['url'],cache,budget=args.budget_mib*1024**2)
            try:
                assert reader.size==spec['bytes']
                with h5py.File(reader,'r') as handle:
                    assert np.array_equal(index(handle['var']),genes)
                    meta=pd.DataFrame({key:column(handle['obs'][key]) for key in ['target_gene','guide_id','batch']})
                    meta['source_row']=np.arange(len(meta));meta['source_cell_id']=index(handle['obs'])
                    assert meta.source_cell_id.is_unique
                    reused=meta.set_index('source_cell_id',drop=False).loc[control_ids].copy()
                    assert (reused.target_gene=='non-targeting').all()
                    assert np.array_equal(reused.batch.astype(str).to_numpy(),controls.batch.astype(str).to_numpy())
                    x=handle['X'];assert x.attrs.get('encoding-type')=='csr_matrix';indptr=x['indptr'][:]
                    chosen=None
                    for per_batch in [16,8,4,2,1]:
                        pert=select_stratified(meta,meta.target_gene.isin(official).to_numpy(),['target_gene','batch'],per_batch)
                        proposed=pd.concat([reused,pert],ignore_index=True).sort_values('source_row').reset_index(drop=True)
                        try:prefetch_v5_sparse(reader,x,indptr,proposed.source_row.to_numpy());chosen=proposed;break
                        except RuntimeError as error:
                            if 'Planned source ranges' not in str(error):raise
                            print('Reduce',split,'target cells/batch',per_batch,str(error),flush=True)
                    if chosen is None:raise RuntimeError('No useful bounded H1 extraction plan')
                    matrix=read_rows(x,indptr,chosen.source_row.to_numpy(),len(genes));meta=chosen
                    meta['public_split']=split
                    reindex=meta.reset_index().set_index('source_cell_id')['index'].loc[control_ids].to_numpy()
                    difference=matrix[reindex]-control_counts
                    if difference.nnz:raise ValueError('Shared control IDs have different original expression counts')
                    sparse.save_npz(directory/'selected_counts.npz',matrix);meta.to_csv(directory/'selected_cells.csv',index=False)
                    record={'source':reader.identity,'target_cells_per_target_batch':per_batch,'selected_cells':len(meta),
                        'target_cells':int((meta.target_gene!='non-targeting').sum()),'targets':int(meta.loc[meta.target_gene!='non-targeting','target_gene'].nunique()),
                        'shared_controls_checked_count_exact':len(control_ids),'shared_controls_retained_in_combined_matrix_once':True,
                        'budget_bytes':reader.budget,'range_bytes_used':reader.used,'downloaded_this_run':reader.downloaded,
                        'counts_sha256':build.sha256(directory/'selected_counts.npz'),'selected_cells_sha256':build.sha256(directory/'selected_cells.csv'),
                        'script_sha256':build.sha256(Path(__file__)),
                        'selection':'Existing 1536 training controls matched by ID and exactly verified, first N target cells per target/batch; largest predeclared N fitting bounded range plan.'}
                    build.write_json(saved,record)
            finally:reader.close()
        assert (matrix[meta.reset_index().set_index('source_cell_id')['index'].loc[control_ids].to_numpy()]-control_counts).nnz==0
        records.append(record);perturb=meta.target_gene!='non-targeting'
        all_meta.append(meta.loc[perturb]);all_counts.append(matrix[np.flatnonzero(perturb.to_numpy())]);del matrix;gc.collect()
    merged=pd.concat(all_meta,ignore_index=True);matrix=sparse.vstack(all_counts,format='csr')
    assert merged.source_cell_id.is_unique and matrix.shape==(len(merged),len(genes))
    sparse.save_npz(OUT/'selected_counts.npz',matrix);merged.to_csv(OUT/'selected_cells.csv',index=False)
    merged.groupby(['public_split','target_gene']).size().rename('selected_cells').to_csv(OUT/'target_coverage.csv')
    report={'original_training_extraction':provenance,'additional_splits':records,'targets':int(merged.loc[merged.target_gene!='non-targeting','target_gene'].nunique()),
        'control_cells':int((merged.target_gene=='non-targeting').sum()),'selected_cells':len(merged),'genes':len(genes),
        'counts_sha256':build.sha256(OUT/'selected_counts.npz'),'selected_cells_sha256':build.sha256(OUT/'selected_cells.csv'),
        'gene_names_sha256':build.sha256(OUT/'gene_names.csv'),'script_sha256':build.sha256(Path(__file__)),
        'scope':'Released public 2025 H1, combined unique cells; source batch splits globally disjoint; no hidden 2026 response data.'}
    build.write_json(OUT/'extraction.json',report);print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
