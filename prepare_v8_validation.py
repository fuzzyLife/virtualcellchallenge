"""Freeze control-defined technical strata and target-independent V8 evaluation cells."""
import json,gc,hashlib
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
from fit_v7_response import raw_inputs,read_cache
from v8_methods import ROOT,A,SEEDS
from v6_methods import moments
import build_submission as build


def split_controls(meta,seed):
    rng=np.random.default_rng(seed);left=[];right=[]
    controls=np.flatnonzero(meta.gene_target.to_numpy()=='Non-Targeting')
    for _,group in meta.iloc[controls].groupby('sample',sort=True):
        ids=rng.permutation(group.index.to_numpy());half=len(ids)//2
        left.extend(ids[:half]);right.extend(ids[half:])
    return np.asarray(left),np.asarray(right)


def global_target_split(target,seed=20268100):
    import hashlib
    return 'calibration' if int(hashlib.sha256(f'{seed}:{target}'.encode()).hexdigest(),16)%2==0 else 'confirmation'


def main():
    A.mkdir(parents=True,exist_ok=True)
    if (A/'validation_plan.json').exists():raise FileExistsError('Frozen V8 plan already exists; preserve it')
    data=read_cache();genes=data['genes'];targets=data['targets'];contexts=[]
    for line in ['H1','HCT116','HEK293T']:
        x,meta,measured,source=raw_inputs(line,genes);left,right=split_controls(meta,SEEDS[0])
        batches=sorted(meta.iloc[left]['sample'].unique())
        profiles=[]
        for batch in batches:
            ids=left[meta.iloc[left]['sample'].to_numpy()==batch]
            block=x[ids];depth=np.asarray(block.sum(axis=1)).ravel()
            profiles.append(np.asarray(block.multiply((1e6/depth)[:,None]).mean(axis=0)).ravel())
        profiles=np.asarray(profiles);columns=np.flatnonzero((profiles.mean(axis=0)>5)&measured)
        columns=columns[np.argsort(np.var(np.log1p(profiles[:,columns]),axis=0))[-1000:]]
        embedding=PCA(n_components=min(6,len(batches)-1),random_state=SEEDS[0]).fit_transform(np.log1p(profiles[:,columns]))
        k=1 if line=='H1' else 2
        labels=np.zeros(len(batches),int) if k==1 else KMeans(n_clusters=k,random_state=SEEDS[0],n_init=10).fit_predict(embedding)
        # Stable cluster labels from their lexicographically smallest batch, independent of outcomes.
        ordering=sorted(range(k),key=lambda c:min(np.asarray(batches)[labels==c]))
        labels=np.array([ordering.index(c) for c in labels]);lookup=dict(zip(batches,labels.tolist()))
        old=json.loads((ROOT/f'experiments/v7_multisource_adaptive/analysis/cell_split_{line}.json').read_text())
        used=set(old['targets'])
        for cluster in range(k):
            context=f'{line}_controlstratum{cluster}'
            selected_batches=[b for b,c in lookup.items() if c==cluster]
            fit=left[meta.iloc[left]['sample'].isin(selected_batches).to_numpy()]
            ref=right[meta.iloc[right]['sample'].isin(selected_batches).to_numpy()]
            if min(len(fit),len(ref))<400:raise ValueError(f'{context}: insufficient unique control cells')
            depths,p,b,_=moments(x[fit]);query=p*1e6
            chosen=[];counts=meta[meta['sample'].isin(selected_batches)].gene_target.value_counts()
            for split in ['calibration','confirmation']:
                available=[t for t in targets if counts.get(t,0)>=60 and query[np.flatnonzero(genes==t)[0]]>5 and global_target_split(t)==split]
                available.sort(key=lambda t:(str(t) in used,hashlib.sha256(f'V8:{t}'.encode()).hexdigest()))
                if len(available)<4:raise ValueError(f'{context}: insufficient {split} targets')
                chosen.extend(available[:4])
            cells=[]
            for seed in SEEDS:
                rng=np.random.default_rng(seed+cluster)
                reference=rng.choice(ref,400,replace=False)
                # Shared template rows across every target and candidate; never read treated batch mixture.
                template=rng.choice(fit,400,replace=False)
                assert not set(template)&set(reference)
                for target in chosen:
                    actual=np.flatnonzero((meta.gene_target.to_numpy()==target)&meta['sample'].isin(selected_batches).to_numpy())
                    treated=np.random.default_rng(seed+int(np.flatnonzero(targets==target)[0])).choice(actual,min(200,len(actual)),replace=False)
                    replicate=np.random.default_rng(seed+31000+int(np.flatnonzero(targets==target)[0])).permutation(actual)
                    half=min(len(replicate)//2,200)
                    cells.append({'target':str(target),'split':global_target_split(target),'seed':seed,
                                  'actual':treated.tolist(),'replicate_left':replicate[:half].tolist(),
                                  'replicate_right':replicate[half:2*half].tolist(),
                                  'template':template.tolist(),'reference':reference.tolist()})
            contexts.append({'context':context,'family':line,'batches':selected_batches,'left':fit.tolist(),'right':ref.tolist(),
                             'cells':cells,'old_v7_targets_reused':sorted(set(chosen)&used),
                             'selected_targets':list(map(str,chosen)),'source':source})
            print('Frozen',context,'controls',len(fit),len(ref),'targets',len(chosen),'old reused',len(set(chosen)&used),flush=True)
        pd.DataFrame({'batch':batches,'control_defined_stratum':labels}).to_csv(A/f'batch_strata_{line}.csv',index=False)
        del x;gc.collect()
    plan={'version':8,'seeds':SEEDS,'contexts':contexts,'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
          'script_sha256':build.sha256(ROOT/'prepare_v8_validation.py'),
          'selection':'Control-only batch PCA/KMeans; globally disjoint target-hash calibration/confirmation; metadata count>=60 and control own CPM>5; prefer targets absent from V7 count tests.',
          'templates':'Same 400 unique control rows for every target within each stratum/seed, independent of actual treated batch proportions; disjoint reference controls.',
          'scope':'Public development data, not pristine blind ground truth; residual within-stratum batch differences remain. 60..200 observed source cells, 400 prediction/reference cells; not official anchors.'}
    build.write_json(A/'validation_plan.json',plan)
    print('Validation plan SHA256',build.sha256(A/'validation_plan.json'),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
