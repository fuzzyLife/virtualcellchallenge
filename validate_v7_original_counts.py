"""Independently compare converted-source counts against original Xaira Batch1 Parquet."""
import argparse,json
import numpy as np,pandas as pd,pyarrow.parquet as pq
from scipy import sparse
import build_submission as build
ROOT=build.ROOT

def check(line):
    directory=ROOT/f'reference/v7_sources/{line}';selected=pd.read_csv(directory/'selected_cells.csv')
    sample=selected[selected['sample']==line+'_Batch1'].head(32)
    assert len(sample)==32
    original=ROOT/f'reference/orion_pilot/{line}_Batch1.parquet';metadata=ROOT/'reference/orion_pilot/gene_metadata.parquet'
    table=pq.read_table(original,filters=[('cell_barcode','in',sample.cell_id.tolist())]).to_pandas().set_index('cell_barcode')
    assert table.index.is_unique and len(table)==len(sample)
    _,genes,_=build.read_spec(ROOT);positions={g:i for i,g in enumerate(genes)}
    var=pq.read_table(metadata).to_pandas();dest={int(r.gene_token_id):positions.get(str(r.gene_name),-1) for r in var.itertuples()}
    chunks={};comparisons=[]
    for r in sample.itertuples():
        source=table.loc[r.cell_id];tokens=np.array(source.gene_token_id,dtype=np.int64);raw=np.array(source.gene_expression)
        assert np.isfinite(raw).all() and (raw>=0).all() and (raw==np.floor(raw)).all()
        assert int(raw.sum())==int(r.total_counts)==int(source.total_counts)
        target=np.array([dest[int(token)] for token in tokens]);keep=target>=0
        expected=sparse.csr_matrix((raw[keep].astype(np.int32),(np.zeros(keep.sum(),dtype=int),target[keep])),shape=(1,len(genes)))
        chunk_number=r.Index//256
        if chunk_number not in chunks:chunks[chunk_number]=sparse.load_npz(directory/'chunks'/f'{chunk_number:05d}.npz')
        actual=chunks[chunk_number][r.Index%256]
        assert (actual!=expected).nnz==0
        comparisons.append({'cell_id':str(r.cell_id),'gene_target':str(r.gene_target),'library':int(raw.sum()),'nnz':int(expected.nnz)})
    record={'passed':True,'line':line,'cells':len(comparisons),'all_18533_gene_counts_match_original_parquet':True,
        'original_parquet_sha256':build.sha256(original),'original_gene_metadata_sha256':build.sha256(metadata),
        'selected_cells_sha256':build.sha256(directory/'selected_cells.csv'),'script_sha256':build.sha256(ROOT/'validate_v7_original_counts.py'),
        'comparisons':comparisons,'scope':'Exact raw count and gene mapping check for 32 distinct original Batch1 cells, not a complete original-release audit.'}
    build.write_json(ROOT/f'analysis/v7_source_inventory/{line}_original_count_replay.json',record)
    print(line,'32 original Batch1 cells: exact complete challenge-axis counts and full libraries verified',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--line',choices=['HCT116','HEK293T','all'],default='all');a=p.parse_args()
    for line in ['HCT116','HEK293T'] if a.line=='all' else [a.line]:check(line)
