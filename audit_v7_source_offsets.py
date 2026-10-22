"""Audit converted Orion expression gaps and repair metadata offsets before extraction."""
import argparse,json
import numpy as np,pandas as pd
from prepare_v7_sources import reader,OUT,REVISION,revision_check
import build_submission as build

def audit(line):
    directory=OUT/line;file=directory/'expression_fragment_audit.json'
    if file.exists():
        cached=json.loads(file.read_text())
        if 'missing_expression_records' in cached:return cached
    revision_check();ds=reader(line,'expression');fs=ds.get_fragments()
    counts=np.array([f.metadata.num_rows for f in fs],dtype=np.int64);starts=np.r_[0,np.cumsum(counts)[:-1]]
    partial=directory/'expression_fragment_audit_partial.json';bounds=[]
    if partial.exists():
        rec=json.loads(partial.read_text());assert rec['counts']==counts.tolist() and rec['start']==starts.tolist();bounds=rec['cell_bounds']
    for k in range(len(bounds),len(fs),128):
        stop=min(k+128,len(fs));probe=np.column_stack([starts[k:stop],starts[k:stop]+counts[k:stop]-1]).ravel()
        bounds.extend(ds.take(probe.tolist(),columns=['cell_integer_id'])['cell_integer_id'].to_numpy().reshape(-1,2).tolist())
        build.write_json(partial,{'counts':counts.tolist(),'start':starts.tolist(),'cell_bounds':bounds})
        print(line,'audited expression fragments',stop,'/',len(fs),flush=True)
    bounds=np.array(bounds);assert np.all(bounds[1:,0]>=bounds[:-1,1])
    gaps=[[int(a+1),int(c-1)] for a,c in zip(bounds[:-1,1],bounds[1:,0]) if c>a+1]
    cells=reader(line,'cells');last=int(cells.count_rows())-1
    if bounds[0,0]>0:gaps.insert(0,[0,int(bounds[0,0])-1])
    if bounds[-1,1]<last:gaps.append([int(bounds[-1,1])+1,last])
    counts_missing=[]
    for lo,hi in gaps:
        rows=cells.to_table(columns=['cell_integer_id','n_genes_by_counts'],filter=f'cell_integer_id >= {lo} AND cell_integer_id <= {hi}').to_pandas()
        assert rows.cell_integer_id.is_unique and len(rows)==hi-lo+1
        counts_missing.append(int(rows.n_genes_by_counts.sum()))
    rec={'revision':REVISION,'counts':counts.tolist(),'start':starts.tolist(),'cell_bounds':bounds.tolist(),
        'missing_cell_ranges':gaps,'missing_expression_records':counts_missing,'expression_rows':int(counts.sum())}
    build.write_json(file,rec);return rec

def repair(line):
    directory=OUT/line;record=audit(line);file=directory/'selected_cells.csv'
    selection=json.loads((directory/'selection.json').read_text())
    if selection.get('expression_offsets_repaired'):
        assert build.sha256(file)==selection['selected_sha256'];return
    original=pd.read_csv(file);old_hash=build.sha256(file);original.to_csv(directory/'selected_cells_unfiltered.csv',index=False)
    keep=np.ones(len(original),bool);offset=np.zeros(len(original),np.int64)
    for (lo,hi),count in zip(record['missing_cell_ranges'],record['missing_expression_records']):
        keep&=~original.cell_integer_id.between(lo,hi).to_numpy();offset[original.cell_integer_id.to_numpy()>hi]+=count
    revised=original.loc[keep].copy();revised['source_cell_start_index']=revised.cell_start_index
    revised['cell_start_index']=revised.cell_start_index-offset[keep];revised=revised.reset_index(drop=True)
    # Read probes distributed across the entire retained cell axis before accepting correction.
    probe=np.unique(np.linspace(0,len(revised)-1,min(128,len(revised))).astype(int));ds=reader(line,'expression')
    ids=ds.take(revised.iloc[probe].cell_start_index.to_numpy(dtype=np.uint64).tolist(),columns=['cell_integer_id'])['cell_integer_id'].to_numpy()
    np.testing.assert_array_equal(ids,revised.iloc[probe].cell_integer_id.to_numpy())
    revised.to_csv(file,index=False);new_hash=build.sha256(file)
    # Reuse a verified chunk only when its exact rows, offsets and source counts did not change.
    reused=0
    for side in sorted((directory/'chunks').glob('*.json')):
        rec=json.loads(side.read_text());a,b=rec['start'],rec['stop'];chunk=side.with_suffix('.npz')
        columns=['cell_integer_id','cell_start_index','n_genes_by_counts','total_counts']
        if rec['selected_sha256']==old_hash and b<=len(revised) and original.iloc[a:b][columns].reset_index(drop=True).equals(revised.iloc[a:b][columns].reset_index(drop=True)) and build.sha256(chunk)==rec['sha256']:
            rec['selected_sha256']=new_hash;rec['offset_repair_reuse_verified']=True;build.write_json(side,rec);reused+=1
    selection.update(selected_sha256=new_hash,selected_cells=len(revised),expression_offsets_repaired=True,
        original_selected_sha256=old_hash,missing_selected_cells=int((~keep).sum()),expression_audit_sha256=build.sha256(directory/'expression_fragment_audit.json'),
        offset_repair='Exclude cells absent from expression fragments; subtract complete missing-cell expression record counts from later global offsets; verify 128 retained cell starts across full axis.',
        reused_verified_chunks=reused,targets=int(revised.loc[revised.gene_target!='Non-Targeting','gene_target'].nunique()))
    build.write_json(directory/'selection.json',selection)
    print(line,'repaired offsets; excluded',int((~keep).sum()),'missing selected cells; reused',reused,'verified chunks',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--line',choices=['HCT116','HEK293T','all'],default='all');args=p.parse_args()
    for line in ['HCT116','HEK293T'] if args.line=='all' else [args.line]:repair(line)
