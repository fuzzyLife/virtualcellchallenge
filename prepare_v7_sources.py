"""Stratified raw CRISPRi cells across the complete pinned Orion release."""
from pathlib import Path
import argparse,json,sys,time,gc,os
from concurrent.futures import ThreadPoolExecutor
import requests,numpy as np,pandas as pd,pyarrow as pa
from scipy import sparse
import pyarrow.parquet as pq
import build_submission as build
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'reference/v7_sources'
REVISION='598aa5442a4ee9a8f1a383a843ab5ce142904d04'
REPO='slaf-project/X-Atlas-Orion'
SEED=20261008
os.environ.setdefault('LANCE_IO_THREADS','8')
def reader(line,axis):
    sys.path.insert(0,str(ROOT/'reference/v7_runtime'))
    import lance
    for attempt in range(12):
        try:
            return lance.dataset(f'hf://datasets/{REPO}/data/{line}/{axis}.lance',storage_options={'hf_revision':REVISION,'hf_enable_resolve_cache':'true','lance_aimd_initial_rate':'100','lance_aimd_max_rate':'300','lance_aimd_burst_capacity':'10'})
        except (ValueError,OSError) as error:
            if attempt==11:raise
            print(line,axis,'retry opening public source',attempt+1,flush=True)
            time.sleep(30)
def revision_check():
    r=requests.get(f'https://huggingface.co/api/datasets/{REPO}',timeout=60);r.raise_for_status()
    if r.json()['sha']!=REVISION:raise ValueError('Atlas release changed; preserve cache')
def balanced_sample(frame,cap,seed):
    rng=np.random.default_rng(seed);f=frame.iloc[rng.permutation(len(frame))].copy()
    f['_within']=f.groupby('sample',sort=False).cumcount()
    return f.sort_values('_within',kind='stable').head(cap).drop(columns='_within')
def metadata(line,cap=400,controls_per_batch=32):
    directory=OUT/line;directory.mkdir(parents=True,exist_ok=True)
    file=directory/'selected_cells.csv'
    if file.exists():
        record=json.loads((directory/'selection.json').read_text())
        assert record['cap']==cap and record['controls_per_batch']==controls_per_batch
        assert build.sha256(file)==record['selected_sha256']
        return pd.read_csv(file)
    revision_check()
    targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.astype(str).to_numpy()
    quoted=','.join("'"+s.replace("'","''")+"'" for s in [*targets,'Non-Targeting'])
    columns=['sample','guide_target','gene_target','total_counts','n_genes_by_counts','pass_guide_filter','cell_id','cell_integer_id','cell_start_index']
    print(line,'reading full-release target/control metadata',flush=True)
    all_cells=reader(line,'cells').to_table(columns=columns,filter=f'gene_target IN ({quoted})').to_pandas()
    assert all_cells.cell_id.is_unique and all_cells.cell_integer_id.is_unique
    all_cells=all_cells[all_cells.pass_guide_filter&(all_cells.total_counts>0)].copy()
    all_cells.groupby('gene_target').size().rename('available_cells').to_csv(directory/'available_target_counts.csv')
    batches=sorted(all_cells['sample'].unique());rng=np.random.default_rng(SEED)
    batches=np.asarray(batches)[rng.permutation(len(batches))]
    partition={str(b):i%2 for i,b in enumerate(batches)}
    all_cells['partition']=all_cells['sample'].map(partition)
    chosen=[]
    for ti,target in enumerate(targets):
        for part in [0,1]:
            group=all_cells[(all_cells.gene_target==target)&(all_cells.partition==part)]
            if len(group):chosen.append(balanced_sample(group,cap//2,SEED+ti*10+part))
    for bi,(batch,group) in enumerate(all_cells[all_cells.gene_target=='Non-Targeting'].groupby('sample',sort=True)):
        chosen.append(balanced_sample(group,controls_per_batch,SEED+3000+bi))
    f=pd.concat(chosen,ignore_index=True).sort_values('cell_integer_id').reset_index(drop=True)
    assert f.cell_id.is_unique
    f.to_csv(file,index=False)
    build.write_json(directory/'selection.json',{'revision':REVISION,'repo':REPO,'cap':cap,'controls_per_batch':controls_per_batch,
        'selected_sha256':build.sha256(file),'available_passing_cells':len(all_cells),'selected_cells':len(f),'batches':len(batches),
        'targets':int(f.loc[f.gene_target!='Non-Targeting','gene_target'].nunique()),'partition_by_batch':partition,
        'selection':'Random within batch; round-robin across batches; cap/2 cells per target per disjoint batch partition; controls maximum 32 per batch.'})
    print(line,'selected',len(f),'cells;',len(batches),'batches',flush=True)
    return f
def extract(line,cap=400,chunk_cells=256,workers=3):
    revision_check();directory=OUT/line;directory.mkdir(parents=True,exist_ok=True);done=directory/'extraction.json'
    canonical_file=ROOT/'reference/orion_pilot/gene_metadata.parquet';mapping_hash=build.sha256(canonical_file)
    if done.exists() and json.loads(done.read_text()).get('canonical_gene_metadata_sha256')==mapping_hash:
        record=json.loads(done.read_text());assert record['revision']==REVISION and record['cap']==cap
        for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_sha256'),('gene_names.csv','genes_sha256')]:
            assert build.sha256(directory/name)==record[key]
        print(line,'verified completed extraction',flush=True);return
    metadata(line,cap)
    from audit_v7_source_offsets import repair
    repair(line)
    selected=pd.read_csv(directory/'selected_cells.csv');selected_hash=build.sha256(directory/'selected_cells.csv')
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.astype(str).to_numpy()
    var=pd.read_csv(directory/'source_genes.csv') if (directory/'source_genes.csv').exists() else reader(line,'genes').to_table(columns=['gene_id','gene_integer_id']).to_pandas()
    var.to_csv(directory/'source_genes.csv',index=False);positions={g:i for i,g in enumerate(genes)}
    canonical=pq.read_table(canonical_file).to_pandas().set_index('gene_token_id')
    assert canonical.index.is_unique and set(var.gene_integer_id)==set(canonical.index)
    # The converted table suffixes duplicate symbols; original integer IDs retain canonical symbols.
    mapping=np.full(int(var.gene_integer_id.max())+1,-1,np.int32)
    for row in var.itertuples():mapping[int(row.gene_integer_id)]=positions.get(str(canonical.loc[int(row.gene_integer_id),'gene_name']),-1)
    measured=np.isin(np.arange(len(genes)),mapping[mapping>=0])
    pd.DataFrame({'gene_name':genes,'measured':measured}).to_csv(directory/'gene_names.csv',index=False)
    ds=reader(line,'expression');chunks=directory/'chunks';chunks.mkdir(exist_ok=True)
    jobs=[(i,start,min(start+chunk_cells,len(selected))) for i,start in enumerate(range(0,len(selected),chunk_cells))]
    def fetch(job):
        number,start,stop=job;file=chunks/f'{number:05d}.npz';sidecar=file.with_suffix('.json')
        if file.exists() and sidecar.exists():
            rec=json.loads(sidecar.read_text())
            if rec['start']==start and rec['stop']==stop and rec['selected_sha256']==selected_hash and build.sha256(file)==rec['sha256'] and rec.get('canonical_gene_metadata_sha256')==mapping_hash:return number
        rows=selected.iloc[start:stop];sizes=rows.n_genes_by_counts.to_numpy(dtype=np.int64)
        indices=np.concatenate([np.arange(int(r.cell_start_index),int(r.cell_start_index)+int(r.n_genes_by_counts),dtype=np.uint64) for r in rows.itertuples()])
        for attempt in range(12):
            try:
                table=ds.take(pa.array(indices),columns=['cell_integer_id','gene_integer_id','value']);break
            except OSError as error:
                print(line,'retry chunk',number,attempt+1,str(error)[:120],flush=True)
                if attempt==11:raise
                time.sleep(min(5*(attempt+1),20))
        ci=table['cell_integer_id'].to_numpy();gi=table['gene_integer_id'].to_numpy()
        raw=table['value'].to_numpy()
        if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any():raise ValueError('Invalid raw counts')
        values=raw.astype(np.int64);np.testing.assert_array_equal(ci,np.repeat(rows.cell_integer_id.to_numpy(),sizes))
        pointers=np.r_[0,np.cumsum(sizes)];libraries=np.add.reduceat(values,pointers[:-1])
        np.testing.assert_array_equal(libraries,rows.total_counts.to_numpy(dtype=np.int64))
        dest=mapping[gi];keep=dest>=0;rr=np.repeat(np.arange(len(rows),dtype=np.int32),sizes)[keep]
        matrix=sparse.csr_matrix((values[keep].astype(np.int32),(rr,dest[keep])),shape=(len(rows),len(genes)))
        matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices();sparse.save_npz(file,matrix)
        build.write_json(sidecar,{'start':start,'stop':stop,'sha256':build.sha256(file),'selected_sha256':selected_hash,
            'full_library_verified':True,'source_records':int(len(values)),'canonical_gene_metadata_sha256':mapping_hash})
        return number
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i,_ in enumerate(pool.map(fetch,jobs)):
            if (i+1)%10==0 or i+1==len(jobs):print(line,'count chunks',i+1,'/',len(jobs),flush=True)
    revision_check();print(line,'assembling verified cell cache',flush=True)
    matrix=sparse.vstack([sparse.load_npz(chunks/f'{job[0]:05d}.npz') for job in jobs],format='csr')
    assert matrix.shape==(len(selected),len(genes));sparse.save_npz(directory/'selected_counts.npz',matrix)
    build.write_json(done,{'revision':REVISION,'repo':REPO,'cap':cap,'cells':len(selected),
        'targets':int(selected.loc[selected.gene_target!='Non-Targeting','gene_target'].nunique()),'batches':int(selected['sample'].nunique()),
        'measured_genes':int(measured.sum()),'counts_sha256':build.sha256(directory/'selected_counts.npz'),'selected_sha256':selected_hash,
        'genes_sha256':build.sha256(directory/'gene_names.csv'),'source_genes_sha256':build.sha256(directory/'source_genes.csv'),
        'script_sha256':build.sha256(Path(__file__)),'offset_audit_sha256':build.sha256(directory/'expression_fragment_audit.json'),'selection_record':json.loads((directory/'selection.json').read_text()),
        'canonical_gene_metadata_sha256':mapping_hash,'canonical_gene_mapping':'Original gene token IDs mapped to original symbols; duplicate source symbols summed on challenge axis, never dropped because of converted suffixes.',
        'normalization':'Full-source library verified before gene-axis projection; missing genes remain masked.',
        'workers':workers,'lance_io_threads':int(os.environ['LANCE_IO_THREADS']),
        'scope':'Stratified target and control sample across full release batches, not every raw cell.'})
    del matrix;gc.collect()
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--line',choices=['HCT116','HEK293T','all'],default='all')
    p.add_argument('--cap',type=int,default=400);p.add_argument('--workers',type=int,default=3);p.add_argument('--metadata-only',action='store_true');p.add_argument('--wait-for-offset-audit',action='store_true')
    a=p.parse_args()
    for line in ['HCT116','HEK293T'] if a.line=='all' else [a.line]:
        if a.metadata_only:metadata(line,a.cap)
        else:
            if a.wait_for_offset_audit:
                while not (OUT/line/'expression_fragment_audit.json').exists():time.sleep(20)
            extract(line,a.cap,workers=a.workers)

