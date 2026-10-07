"""Fit a small embedding-assisted model of public CRISPRi response contrasts.

Run: .\\.competition_env\\Scripts\\python.exe train_embedding_response.py
All source-lineage and target holdouts are evaluated. No leaderboard upload.
"""
from pathlib import Path
import json, gc
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from threadpoolctl import threadpool_limits
import build_submission as build
from build_cd4_model import eligible, load_subset
from fit_context_ensemble import target_folds
from fit_lineage_submission import bh_valid
from embedding_response_methods import VARIANTS, RIDGE, kernels, predict_ridge

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'analysis/embedding_response_pilot'
CONTEXTS = np.array(['K562', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr', 'H1'])
FAMILIES = np.array(['K562', 'CD4', 'CD4', 'CD4', 'H1'])


def embed_file(name, genes):
    path = ROOT/'cross_modal_foundation_model/phase1/data/phase1/benchmarks/depmap/gene-embs'/name
    with h5py.File(path) as f:
        symbols = f['obs/gene'].asstr()[:]
        raw = f['X'][:].astype(float)
    if len(set(symbols)) != len(symbols) or not np.isfinite(raw).all():
        raise ValueError('Invalid embedding axis/values')
    mapping = {g: i for i, g in enumerate(symbols)}
    vectors = np.zeros((len(genes), raw.shape[1]))
    available = np.array([g in mapping for g in genes])
    vectors[available] = raw[[mapping[g] for g in genes[available]]]
    return (vectors, available), path


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    _, genes, targets = build.read_spec(ROOT)
    gp = {g: i for i, g in enumerate(genes)}
    tp = {t: i for i, t in enumerate(targets)}
    outputs = ~np.isin(genes, targets)
    profiles, hashes = [], {}
    for context in CONTEXTS:
        path = ROOT/f'reference/control_profiles/{context}.csv'
        record = json.loads(path.with_suffix('.json').read_text())
        assert build.sha256(path) == record['profile_sha256']
        frame = pd.read_csv(path).set_index('gene_name')
        assert frame.index.is_unique
        profiles.append(frame.pooled_cpm.reindex(genes).to_numpy())
        hashes[str(path.relative_to(ROOT))] = build.sha256(path)
    profiles = np.asarray(profiles)
    values, truths, errors, qvalues, masks, row_records = [], [], [], [], [], []
    def append(context, names, sg, y, se, q, n, estimate):
        si = np.array([i for i, g in enumerate(sg) if g in gp]);di=np.array([gp[sg[i]] for i in si])
        # Coverage uses gene identity and controls only. Label availability is
        # masked per row, so a held-out target cannot change the training mask.
        covered = outputs[di].copy()
        if context == 4:
            covered &= np.isfinite(profiles[context,di]) & (profiles[context,di]>5)
        for i, name in enumerate(names):
            if name not in tp:
                continue
            usable = covered & np.isfinite(y[i,si]) & np.isfinite(se[i,si]) & (se[i,si]>0)
            v = np.zeros(len(genes),bool);v[di[usable]]=True
            raw=np.zeros(len(genes));raw[di]=np.nan_to_num(np.clip(y[i,si],-3,3))
            est=np.zeros(len(genes));est[di]=np.nan_to_num(np.clip(estimate[i,si],-3,3));est[~v]=0
            err=np.full(len(genes),1e6);err[di]=np.nan_to_num(se[i,si],nan=1e6)
            qq=np.ones(len(genes));qq[di]=np.nan_to_num(q[i,si],nan=1.)
            values.append(est);truths.append(raw);errors.append(err);qvalues.append(qq);masks.append(v)
            row_records.append({'row':len(row_records),'context':CONTEXTS[context], 'context_index':context,
                'family':FAMILIES[context], 'target_gene':name,'target_index':tp[name],
                'source_cells':int(n[i]),'modeled_genes':int(v.sum())})
    # K562: observed pseudobulk contrasts; keep original effective-guide selection.
    path=ROOT/'reference/K562_gwps_raw_bulk_01.h5ad'
    old=json.loads((ROOT/'experiments/v4_lineage_response/analysis/k562_source/reference_model.json').read_text())
    assert build.sha256(path)==old['source_file_sha256'];hashes[str(path.relative_to(ROOT))]=old['source_file_sha256']
    with h5py.File(path) as f:
        labels=f['obs/gene_transcript'].asstr()[:];nc=f['obs/num_cells_filtered'][:].astype(float)
        fold=f['obs/fold_expr'][:];core=f['obs/core_control'][:].astype(bool)
        col=f['var/gene_name'];symbols=f[col.attrs['categories']].asstr()[:][col[:]];matrix=f['X'][:].astype(float)
    sg=np.array(sorted(set(symbols)&set(genes)));sp={g:i for i,g in enumerate(sg)}
    cols=np.flatnonzero(np.isin(symbols,sg))
    mapping=sparse.csr_matrix((np.ones(len(cols)),(cols,[sp[symbols[i]] for i in cols])),shape=(len(symbols),len(sg)))
    matrix=mapping.T.dot(matrix.T).T
    names=np.array([s.split('_',2)[1] if len(s.split('_',2))==3 and s.split('_',1)[0].isdigit() else '' for s in labels])
    ctrl=(names=='non-targeting')&core&(nc>=30)
    base=np.average(matrix[ctrl],axis=0,weights=nc[ctrl]);bcpm=base*1e6/base.sum()
    yr, er, nr, nt=[],[],[],[]
    for target in targets:
        rows=np.flatnonzero((names==target)&(nc>=30)&np.isfinite(fold)&(fold<=.5))
        if not len(rows):continue
        mean=np.average(matrix[rows],axis=0,weights=nc[rows]);cpm=mean*1e6/mean.sum();n=nc[rows].sum()
        y=np.clip(np.log2((cpm+1)/(bcpm+1)),-3,3)
        se=np.sqrt(1/np.maximum(n*mean,1)+1/np.maximum(nc[ctrl].sum()*base,1))/np.log(2)
        yr.append(y);er.append(se);nr.append(n);nt.append(target)
    yr,er=np.asarray(yr),np.asarray(er)
    kk=bh_valid(2*stats.norm.sf(np.abs(yr)/er),np.broadcast_to(bcpm>5,yr.shape))
    append(0,nt,sg,yr,er,kk,nr,np.sign(yr)*np.maximum(np.abs(yr)-1.96*er,0))
    del matrix,mapping;gc.collect()
    # CD4 uses merged donor summaries, retaining the three culture conditions.
    path=ROOT/'reference/cd4_probe/merged.npz'
    v4meta=json.loads((ROOT/'experiments/v4_lineage_response/analysis/reference_model.json').read_text())
    pinned = {key.replace('\\', '/'): value for key, value in v4meta['source_data_sha256'].items()}
    for name in ['reference/cd4_probe/merged.npz','reference/cd4_probe/merged.csv']:
        assert build.sha256(ROOT/name)==pinned[name]
        hashes[name]=pinned[name]
    meta,data=load_subset(path.parent,'merged');good=eligible(meta)
    for ci,condition in enumerate(['Rest','Stim8hr','Stim48hr'],1):
        rows=np.flatnonzero(good&(meta.culture_condition.to_numpy()==condition))
        y=data['log_fc'][rows].astype(float);se=data['lfcSE'][rows].astype(float)
        append(ci,meta.iloc[rows].target_contrast_gene_name.to_numpy(),data['genes'],y,se,
               data['adj_p_value'][rows],meta.iloc[rows].n_cells_target.to_numpy(), y*.04/(.04+se**2))
    # H1 full contrast summaries are labels only when H1 is a training lineage.
    path=ROOT/'experiments/v5_atlas_calibrated/analysis/h1_contrasts.npz'
    cm=json.loads(path.with_suffix('.json').read_text());assert build.sha256(path)==cm['arrays_sha256']
    hashes[str(path.relative_to(ROOT))]=cm['arrays_sha256']
    with np.load(path,allow_pickle=False) as z:
        y=z['2_values'];se=z['2_errors']
        append(4,z['targets'],z['genes'],y,se,z['2_q'],z['2_n'],.7*np.sign(y)*np.maximum(np.abs(y)-se,0))
    emb={}
    for name,file in [('tahoe','tx1_70m-mean-lt5gt70-bin.h5ad'),('pca','pcloading-lt5gt70-bin.h5ad')]:
        emb[name],path=embed_file(file,genes);hashes[str(path.relative_to(ROOT))]=build.sha256(path)
    records=pd.DataFrame(row_records)
    assert not records.duplicated(['context','target_gene']).any(), 'Repeated contrasts require explicit aggregation'
    records.to_csv(OUT/'source_rows.csv',index=False)
    arrays=dict(genes=genes.astype(str),targets=targets.astype(str),profiles=profiles,
                values=np.asarray(values,np.float32),truth=np.asarray(truths,np.float32),errors=np.asarray(errors,np.float32),
                q=np.asarray(qvalues,np.float32),valid=np.asarray(masks),
                row_targets=records.target_index.to_numpy(),row_contexts=records.context_index.to_numpy(),
                row_families=records.family.to_numpy(dtype=str))
    for name,(vectors,available) in emb.items():arrays[name+'_vectors']=vectors;arrays[name+'_available']=available
    np.savez_compressed(OUT/'source_cache.npz',**arrays)
    build.write_json(OUT/'source_provenance.json',{'sources':hashes,'cache_sha256':build.sha256(OUT/'source_cache.npz'),
        'script_sha256':build.sha256(Path(__file__)),'row_counts':records.groupby('context').size().to_dict(),
        'label_units':'clipped log2 response contrasts; fixed source-specific noise shrinkage for regression labels',
        'K562_q_scope':'Poisson-noise normal approximation, not actual cell-level rank DE',
        'H1_q_scope':'Existing source-cell Wilcoxon/BH calls',
        'CD4_q_scope':'Published adjusted p-values from source summaries',
        'mask_scope':'Drop all panel output genes; H1 coverage uses non-targeting control CPM >5. Per-row unavailable responses/errors are masked, never zero-imputed into the loss.'})
    print('Source cache:',records.groupby('context').size().to_dict(),flush=True)
    return arrays


def score(prediction, data, rows, supported, variant, regime, fold):
    result=[];genes,targets=data['genes'],data['targets']
    for i,row in enumerate(rows):
        context=data['row_contexts'][row];target=targets[data['row_targets'][row]]
        valid=data['valid'][row]&supported&np.isfinite(data['profiles'][context])&(data['profiles'][context]>5)
        if not valid.any():continue
        truth=data['truth'][row];real=valid&(data['q'][row]<.05)
        p=2*stats.norm.sf(np.abs(prediction[i])/np.maximum(data['errors'][row],1e-12))
        called=valid&(bh_valid(p[None,:],valid[None,:])[0]<.05)
        correct=called&(np.sign(prediction[i])==np.sign(truth))
        union=called|real;budget=max(called.sum(),real.sum())
        result.append({'regime':regime,'fold':fold,'variant':variant,'family':str(data['row_families'][row]),
            'context':str(CONTEXTS[context]),'target_gene':str(target),'target_index':int(data['row_targets'][row]),
            'mse':float(np.mean((prediction[i,valid]-truth[valid])**2)),
            'zero_mse':float(np.mean(truth[valid]**2)),
            'fidelity_proxy':float(correct.sum()/budget) if budget else 0.,
            'jaccard_proxy':float((called&real).sum()/union.sum()) if union.any() else 1.,
            'sign_accuracy_reference_DE':float((np.sign(prediction[i,real])==np.sign(truth[real])).mean()) if real.any() else None,
            'predicted_calls':int(called.sum()),'reference_calls':int(real.sum()),'genes_evaluated':int(valid.sum())})
    return result


def fit():
    path=OUT/'source_cache.npz';record=json.loads((OUT/'source_provenance.json').read_text())
    assert build.sha256(path)==record['cache_sha256']
    with np.load(path,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    embeddings={k:(data[k+'_vectors'][[np.flatnonzero(data['genes']==t)[0] for t in data['targets']]],
                   data[k+'_available'][[np.flatnonzero(data['genes']==t)[0] for t in data['targets']]]) for k in ['pca','tahoe']}
    gene_vectors=(data['tahoe_vectors'],data['tahoe_available'])
    folds=target_folds(data['targets'],3);scores=[];split_records=[];pred_cache={}
    splits=[]
    for family in ['K562','CD4','H1']:
        test=data['row_families']==family;splits.append(('lineage_holdout',family,~test,test))
    for fold in range(3):
        test=folds[data['row_targets']]==fold;splits.append(('target_holdout',str(fold),~test,test))
    for regime,fold,train,test in splits:
        tr,te=np.flatnonzero(train),np.flatnonzero(test)
        if regime=='lineage_holdout':assert not set(data['row_families'][tr])&set(data['row_families'][te])
        else:assert not set(data['row_targets'][tr])&set(data['row_targets'][te])
        # Fixed evaluation universe: available in at least one TRAINING source,
        # independent of representation. Coverage is evaluated separately.
        supported=data['valid'][tr].any(axis=0)
        split_records.append({'regime':regime,'fold':fold,'train_rows':tr.tolist(),'test_rows':te.tolist(),
            'train_targets':sorted(set(data['targets'][data['row_targets'][tr]].tolist())),
            'test_targets':sorted(set(data['targets'][data['row_targets'][te]].tolist())),
            'train_families':sorted(set(data['row_families'][tr].tolist())), 'test_families':sorted(set(data['row_families'][te].tolist()))})
        for variant in VARIANTS:
            k,cross=kernels(data['row_targets'][tr],data['row_contexts'][tr],data['row_targets'][te],
                data['row_contexts'][te],data['profiles'],embeddings,variant,gene_vectors)
            pred,learned=predict_ridge(k,cross,data['values'][tr],data['valid'][tr],data['row_families'][tr])
            assert np.all(learned==supported)
            scores.extend(score(pred,data,te,supported,variant,regime,fold))
            if regime=='lineage_holdout' and fold=='H1':
                pred_cache[variant]=pred;pred_cache[variant+'_rows']=te;pred_cache[variant+'_supported']=supported
            pd.DataFrame(scores).to_csv(OUT/'contrast_scores.csv',index=False)
            print(regime,fold,variant,'rows',len(te),'done',flush=True)
    frame=pd.DataFrame(scores)
    summary=frame.groupby(['regime','family','variant'])[['mse','zero_mse','fidelity_proxy','jaccard_proxy',
        'sign_accuracy_reference_DE','predicted_calls','reference_calls','genes_evaluated']].mean().reset_index()
    summary['mse_ratio_to_zero']=summary.mse/summary.zero_mse
    summary.to_csv(OUT/'contrast_summary.csv',index=False)
    np.savez_compressed(OUT/'h1_lineage_predictions.npz',**pred_cache)
    build.write_json(OUT/'split_manifest.json',split_records)
    build.write_json(OUT/'fit_manifest.json',{'ridge':RIDGE,'variants':list(VARIANTS),
        'hyperparameter_policy':'Fixed before evaluation; no outcome-based hyperparameter search',
        'source_cache_sha256':record['cache_sha256'],'prediction_cache_sha256':build.sha256(OUT/'h1_lineage_predictions.npz'),
        'script_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/'embedding_response_methods.py'),
        'limitations':['Small, previously inspected public-source panel; not hidden A/B/C validation.',
            'Summary fidelity/jaccard are normal-call diagnostics, not leaderboard metrics.',
            'Only non-panel downstream genes are fitted. All 300 target output columns use production fallback.',
            'The precomputed vectors are RNA-only Tahoe, not proteomics-pretrained Trove.',
            'Full public response summaries, not independent donors, train within each source family.']})
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--prepare-only',action='store_true');parser.add_argument('--reuse-cache',action='store_true');args=parser.parse_args()
    with threadpool_limits(limits=2):
        if not args.reuse_cache:prepare()
        if not args.prepare_only:fit()
