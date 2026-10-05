"""Bounded anonymous metadata coverage probe of public 2025 H1 splits and promoter annotation."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import re
import h5py
import numpy as np
import pandas as pd
import requests
from prepare_reference_controls import RangeReader,column,index
from build_submission import sha256,write_json
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/upstream_review'
CACHE=ROOT/'reference/upstream_review'
ASSETS=json.loads((ROOT/'Virtual-Cell-Challenge-2026/sources.json').read_text())


def coverage():
    official=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy();official_genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    records=[];sources=[];control_identifiers={}
    for split in ['train','validation','test']:
        spec=ASSETS['h1-'+split]
        cache=ROOT/'reference/control_profiles/ranges_H1' if split=='train' else CACHE/('ranges_H1_'+split)
        reader=RangeReader(spec['url'],cache,budget=128*1024**2)
        try:
            assert reader.size==spec['bytes']
            with h5py.File(reader,'r') as f:
                labels=column(f['obs/target_gene']);genes=index(f['var'])
                identifiers=index(f['obs']);assert len(identifiers)==len(labels)
                control_identifiers[split]=set(identifiers[labels=='non-targeting'].astype(str))
                unique,counts=np.unique(labels,return_counts=True);counts=dict(zip(unique,counts))
                for target in official:
                    records.append({'split':split,'target_gene':target,'available_cells':int(counts.get(target,0)),
                                    'gene_in_source':bool(target in set(genes))})
                sources.append({'split':split,'source':reader.identity,'cells':len(labels),
                    'targets':len(set(labels)-{'non-targeting'}),'control_cells':int(counts.get('non-targeting',0)),
                    'overlapping_2026_targets':len((set(labels)-{'non-targeting'})&set(official)),
                    'genes':len(genes),'shared_2026_genes':len(set(genes)&set(official_genes)),
                    'range_bytes_used':reader.used,'downloaded_this_run':reader.downloaded,
                    'scope':'Observation labels and gene axis only; expression matrix never accessed; counts are available cells, not QC-approved usable cells.'})
            print('H1',split,'metadata overlap',sources[-1]['overlapping_2026_targets'],'targets;',reader.used,'range bytes',flush=True)
        finally:reader.close()
    frame=pd.DataFrame(records);frame.to_csv(OUT/'h1_public_split_coverage.csv',index=False)
    combined=frame.pivot(index='target_gene',columns='split',values='available_cells').fillna(0).astype(int)
    combined['all_splits_cells']=combined.sum(axis=1);combined['new_beyond_training']=(combined.train==0)&(combined.all_splits_cells>0)
    combined.to_csv(OUT/'h1_public_target_coverage.csv')
    write_json(OUT/'h1_coverage_probe.json',{'sources':sources,'training_overlap':int((combined.train>0).sum()),
        'all_splits_overlap':int((combined.all_splits_cells>0).sum()),'additional_targets':combined.index[combined.new_beyond_training].tolist(),
        'total_downloaded_this_run':sum(s['downloaded_this_run'] for s in sources),
        'shared_control_cell_identifiers':{split:len(control_identifiers['train']&control_identifiers[split]) for split in ['validation','test']},
        'control_identifier_scope':'Shared identifiers flag potentially reused controls; count equality was not checked. Do not treat public splits as independent control replicates.',
        'script_sha256':sha256(Path(__file__)),'no_billing_project_or_authentication_used':True})


def promoter_pairs():
    spec=ASSETS['gencode47'];path=CACHE/spec['filename'];path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        pending=path.with_suffix(path.suffix+'.partial');size=0;digest=hashlib.sha256()
        with requests.get(spec['url'],stream=True,timeout=(15,90)) as response:
            response.raise_for_status()
            with pending.open('wb') as handle:
                for block in response.iter_content(1024**2):
                    size+=len(block)
                    if size>spec['bytes']:raise ValueError('Annotation exceeds pinned size')
                    digest.update(block);handle.write(block)
        assert size==spec['bytes'] and digest.hexdigest()==spec['sha256'];pending.replace(path)
    assert path.stat().st_size==spec['bytes'] and sha256(path)==spec['sha256']
    rows=[]
    with gzip.open(path,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if fields[2]!='gene':continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            if 'gene_name' in attrs:
                rows.append({'gene':attrs['gene_name'],'chromosome':fields[0],'strand':fields[6],
                             'tss':int(fields[3] if fields[6]=='+' else fields[4])})
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy();targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    frame=pd.DataFrame(rows);unique=frame[~frame.gene.duplicated(keep=False)].set_index('gene');measured=unique.loc[unique.index.intersection(genes)]
    pairs=[]
    for target in targets:
        if target not in unique.index:continue
        a=unique.loc[target];near=measured[(measured.chromosome==a.chromosome)&((measured.tss-a.tss).abs()<=5000)]
        for name,b in near.iterrows():
            if target!=name:pairs.append({'target':target,'neighbor':name,'distance':abs(int(a.tss)-int(b.tss)),
                'chromosome':a.chromosome,'target_tss':int(a.tss),'neighbor_tss':int(b.tss),
                'target_strand':a.strand,'neighbor_strand':b.strand})
    pairs=pd.DataFrame(pairs);pairs.to_csv(OUT/'promoter_neighbor_pairs.csv',index=False)
    meta=pd.read_csv(ROOT/'reference/cd4_probe/merged.csv')
    with np.load(ROOT/'reference/cd4_probe/merged.npz',allow_pickle=False) as data:
        names=data['genes'];lfc=data['log_fc'];q=data['adj_p_value'];positions={str(g):i for i,g in enumerate(names)}
    observed=[]
    for pair in pairs.itertuples():
        if pair.neighbor not in positions:continue
        gi=positions[pair.neighbor]
        for i in np.flatnonzero(meta.target_contrast_gene_name.to_numpy()==pair.target):
            row=meta.iloc[i]
            if row.n_cells_target<100 or not row.ontarget_significant or row.single_guide_estimate or row.distal_offtarget_flag or row.low_target_gex:continue
            observed.append({'target':pair.target,'neighbor':pair.neighbor,'distance':pair.distance,'condition':row.culture_condition,
                'log2fc':float(lfc[i,gi]),'q_value':float(q[i,gi]),'publisher_neighbor_KD_flag':bool(row.neighboring_gene_KD)})
    observed=pd.DataFrame(observed);observed.to_csv(OUT/'cd4_promoter_neighbor_observations.csv',index=False)
    finite=np.isfinite(observed.log2fc)&np.isfinite(observed.q_value);sig=finite&(observed.q_value<.05)
    write_json(OUT/'promoter_probe.json',{'annotation':spec,'pairs':len(pairs),'targets_with_neighbors':pairs.target.nunique(),
        'cd4_observations':len(observed),'significant_negative':int((sig&(observed.log2fc<0)).sum()),
        'significant_positive':int((sig&(observed.log2fc>0)).sum()),
        'scope':'Gene-boundary TSS, not active promoter or guide position; repeated target-condition observations; descriptive source evidence, not challenge validation or permission to impose universal knockdown.',
        'script_sha256':sha256(Path(__file__))})
    print('Promoter pairs',len(pairs),'over',pairs.target.nunique(),'targets; source significant negative/positive',
          int((sig&(observed.log2fc<0)).sum()),int((sig&(observed.log2fc>0)).sum()),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--metadata-only',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);coverage()
    if not args.metadata_only:promoter_pairs()
