"""Bounded, pinned extraction of public 2025 H1 responses for overlapping targets."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from prepare_reference_controls import ROOT, RangeReader, SOURCES, column, index, select_stratified, prefetch_sparse
from build_submission import sha256, write_json

OUT=ROOT/'reference/h1_response'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/'selected_counts.npz'
    if path.exists() and (OUT/'extraction.json').exists():
        report=json.loads((OUT/'extraction.json').read_text())
        if sha256(path)!=report['counts_sha256']:
            raise ValueError('H1 extracted count checksum mismatch')
        print('Verified existing H1 response subset',flush=True); return
    reader=RangeReader(SOURCES['H1'],ROOT/'reference/control_profiles/ranges_H1',budget=2*1024**3)
    if reader.size!=15482497461:
        raise ValueError('Pinned H1 file size changed')
    with h5py.File(reader,'r') as f:
        meta=pd.DataFrame({key:column(f['obs'][key]) for key in ['target_gene','guide_id','batch']})
        meta['source_row']=np.arange(len(meta))
        genes=index(f['var'])
        targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
        mask=meta.target_gene.isin(targets).to_numpy()
        controls=select_stratified(meta,meta.target_gene.to_numpy()=='non-targeting',['batch'],32)
        x=f['X']; indptr=x['indptr'][:]
        chosen=None
        for per_batch in [8,4,2,1]:
            perturb=select_stratified(meta,mask,['target_gene','batch'],per_batch)
            proposed=pd.concat([controls,perturb]).sort_values('source_row')
            try:
                prefetch_sparse(reader,x,indptr,proposed.source_row.to_numpy())
                chosen=proposed; break
            except RuntimeError as error:
                if 'Planned source ranges' not in str(error):
                    raise
                print('Reducing H1 target sampling:',per_batch,str(error),flush=True)
        if chosen is None:
            raise RuntimeError('No H1 response sampling plan fits the two-GiB byte budget')
        data=[]; indices=[]; pointer=[0]
        for count,row in enumerate(chosen.source_row):
            lo,hi=int(indptr[row]),int(indptr[row+1])
            values=x['data'][lo:hi].astype(np.int32); cols=x['indices'][lo:hi].astype(np.int32)
            if (values<0).any() or values.sum(dtype=np.int64)<=0:
                raise ValueError('Invalid H1 counts')
            data.append(values); indices.append(cols); pointer.append(pointer[-1]+len(values))
            if (count+1)%256==0:
                print('H1 rows extracted',count+1,'/',len(chosen),flush=True)
        matrix=sparse.csr_matrix((np.concatenate(data),np.concatenate(indices),np.array(pointer,dtype=np.int64)),shape=(len(chosen),len(genes)))
        matrix.sum_duplicates(); matrix.eliminate_zeros()
        sparse.save_npz(path,matrix)
        chosen.to_csv(OUT/'selected_cells.csv',index=False)
        pd.DataFrame({'gene_name':genes}).to_csv(OUT/'gene_names.csv',index=False)
        coverage=chosen.groupby('target_gene').size().rename('selected_cells').to_frame()
        coverage['available_cells']=meta.groupby('target_gene').size()
        coverage.to_csv(OUT/'target_coverage.csv')
        report={'source':reader.identity,'public_release':'2025 training only','selection':'First source-ordered 32 controls per batch; first N target cells per target/batch, adaptive bounded plan',
                'target_cells_per_batch':per_batch,'targets':int(chosen.loc[chosen.target_gene!='non-targeting','target_gene'].nunique()),'selected_cells':len(chosen),
                'control_cells':len(controls),'genes':len(genes),'range_bytes_used':reader.used,'downloaded_this_run':reader.downloaded,
                'byte_budget':reader.budget,'counts_sha256':sha256(path),'selected_cells_sha256':sha256(OUT/'selected_cells.csv'),
                'gene_names_sha256':sha256(OUT/'gene_names.csv'),'script_sha256':sha256(Path(__file__))}
        write_json(OUT/'extraction.json',report)
        print(json.dumps(report,indent=2),flush=True)
    reader.close()


if __name__=='__main__':
    main()
