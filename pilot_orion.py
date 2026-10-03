"""Bounded pilot of two pinned X-Atlas/Orion shards; never downloads the full atlas.

The pilot measures control-profile similarity and usable target counts, not
perturbation prediction performance. Its sparse per-target samples are not enough
for a full source model. Original challenge artifacts are never overwritten.
"""
from pathlib import Path
import argparse
import gc
import json
import urllib.request
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import build_submission as build

REVISION='53a5bc98d49247bcf967500292575c3d3602de31'
FILES={
'metadata/gene_metadata.parquet':(881915,'33900a0dbaeafb84601c7c8ba7c6642c440774fc3f49613dab11cf1f2bf0fdbf'),
'data/HCT116_Batch1.parquet':(327620255,'a7a6918ecda30343e3f9262bcfb654bf077170d2e99c7bcd0c8b9f963764f13b'),
'data/HEK293T_Batch1.parquet':(370263016,'bc97ddfeef8df72fbd7aa03cf8e5fd1a6293e0d467735e5e59ce907df49af8d1')}


def download(path,name,size,digest):
    if path.exists() and path.stat().st_size==size and build.sha256(path)==digest:
        return
    partial=path.with_suffix('.parquet.partial')
    url=f'https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/{REVISION}/{name}'
    request=urllib.request.Request(url)
    build.log(f'Pilot download: {name} ({size/1e6:.1f} MB)')
    received=0
    with urllib.request.urlopen(request,timeout=60) as response,partial.open('wb') as output:
        while received<size:
            chunk=response.read(min(4*1024*1024,size-received))
            if not chunk:raise RuntimeError('Truncated pilot download')
            output.write(chunk);received+=len(chunk)
        if response.read(1):raise RuntimeError('Pilot file exceeds pinned size')
    if build.sha256(partial)!=digest:raise RuntimeError('Pinned SHA-256 mismatch')
    partial.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();root=args.root.resolve()
    cache=root/'reference/orion_pilot';cache.mkdir(exist_ok=True)
    output=root/'analysis/model_audit';output.mkdir(exist_ok=True)
    _,genes,targets=build.read_spec(root)
    for name,(size,digest) in FILES.items():
        download(cache/Path(name).name,name,size,digest)
    gene_meta=pq.read_table(cache/'gene_metadata.parquet').to_pandas()
    required={'gene_name','gene_token_id'}
    if not required.issubset(gene_meta):raise ValueError(f'Unexpected gene metadata schema: {gene_meta.columns.tolist()}')
    positions={gene:i for i,gene in enumerate(genes)}
    token_to_dest={int(row.gene_token_id):positions[row.gene_name] for row in gene_meta.itertuples() if row.gene_name in positions}
    token_count=int(gene_meta['gene_token_id'].max())+1
    destination=np.full(token_count,-1,dtype=np.int64)
    for token,pos in token_to_dest.items():destination[token]=pos
    root_profiles={c:pd.read_csv(root/f'analysis/gene_statistics_{c}.csv')['mean_cpm'].to_numpy() for c in ['A','B','C']}
    summary={'dataset':'Xaira-Therapeutics/X-Atlas-Orion','revision':REVISION,
             'bytes_downloaded_or_reused':sum(size for size,_ in FILES.values()),
             'file_sha256':{Path(name).name:digest for name,(_,digest) in FILES.items()},
             'gene_metadata_columns':gene_meta.columns.tolist(),'shared_gene_symbols':len(set(token_to_dest.values())),
             'license':'CC BY-NC-SA 4.0','source':'https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion',
             'scope':'One representative sequencing batch from each cell line, not the complete atlas.',
             'limitations':['Guide-library coverage does not imply adequate response coverage in this pilot.',
                            'Shared-symbol comparison does not resolve every historical gene alias.',
                            'Basal similarity alone does not validate transfer of perturbation responses.'], 'cell_lines':{}}
    all_target_counts=[]
    for line in ['HCT116','HEK293T']:
        path=cache/f'{line}_Batch1.parquet'
        build.log(f'Reading {line} pilot counts')
        table=pq.read_table(path,columns=['gene_token_id','gene_expression','gene_target','total_counts','pass_guide_filter'])
        labels=np.array(table.column('gene_target').to_pylist())
        passed=np.array(table.column('pass_guide_filter').to_pylist())>0
        libraries=np.array(table.column('total_counts').to_pylist(),dtype=float)
        control=np.array([str(label).lower().startswith('non-targeting') or str(label).lower() in ['ntc','control'] for label in labels])&passed&(libraries>0)
        if not control.any():raise ValueError(f'No recognizable controls: most common labels {pd.Series(labels).value_counts().head(5).to_dict()}')
        profile=np.zeros(len(genes)); maximum_library_discrepancy=0.
        for row in np.flatnonzero(control):
            tokens=np.array(table.column('gene_token_id')[int(row)].as_py(),dtype=int)
            values=np.array(table.column('gene_expression')[int(row)].as_py(),dtype=float)
            assert len(tokens)==len(values) and np.isfinite(values).all() and (values>=0).all()
            assert (values==np.floor(values)).all(), 'Pilot expressions are not raw whole counts'
            maximum_library_discrepancy=max(maximum_library_discrepancy,float(abs(values.sum()-libraries[row])))
            dest=destination[tokens];keep=dest>=0
            profile+=np.bincount(dest[keep],weights=values[keep]*1e6/libraries[row],minlength=len(genes))
        profile/=control.sum()
        shared=destination[destination>=0]
        shared=np.unique(shared)
        correlations={}
        for context,baseline in root_profiles.items():
            a=profile[shared]*1e6/profile[shared].sum();b=baseline[shared]*1e6/baseline[shared].sum()
            correlations[context]=float(np.corrcoef(np.log1p(a),np.log1p(b))[0,1])
        eligible=passed&(libraries>0)
        counts=pd.Series(labels[eligible]).value_counts()
        target_counts=[{'cell_line':line,'target_gene':target,'pilot_cells':int(counts.get(target,0))} for target in targets]
        all_target_counts.extend(target_counts)
        summary['cell_lines'][line]={'pilot_cells':len(labels),'control_cells':int(control.sum()),
                                    'targets_with_any_cells':sum(row['pilot_cells']>0 for row in target_counts),
                                    'targets_with_at_least_30_cells':sum(row['pilot_cells']>=30 for row in target_counts),
                                    'context_log_cpm_correlations':correlations,
                                    'maximum_control_count_sum_discrepancy':maximum_library_discrepancy}
        pd.DataFrame({'gene_name':genes,'pilot_control_mean_cpm':profile}).to_csv(output/f'orion_{line}_pilot_control_profile.csv',index=False)
        build.log(f'{line}: {summary["cell_lines"][line]}')
        del table;gc.collect()
    pd.DataFrame(all_target_counts).to_csv(output/'orion_pilot_target_counts.csv',index=False)
    build.write_json(output/'orion_pilot_summary.json',summary)
    build.log('Bounded Orion pilot complete; no new prediction generated or submitted')


if __name__=='__main__':main()
