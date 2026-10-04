"""Fit a competition-only shared-gene rank space and project external controls."""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cluster_full_rank_matrix import rank_unit_rows, cluster
from extract_external_mapping_controls import CACHE, load_profiles
from analyze_target_rank_correlations import ROOT, REFERENCES, digest, table

OUT = ROOT / 'analysis/external_rank_mapping'
CONTEXTS = ['A', 'B', 'C']
COLORS = ['#4477aa', '#ee7733', '#228833']


def frozen_projection(raw, model):
    """No refit: apply the same ranks, centroids and PCA to a new measured row."""
    z = rank_unit_rows(raw)
    similarities = z @ model['centroids'].T
    labels = similarities.argmax(axis=1)
    ordered = np.sort(similarities, axis=1)
    pcs = (z - model['pca_mean']) @ model['pca_components'].T
    result = pd.DataFrame({'nearest_context': np.asarray(CONTEXTS)[labels],
                           'centroid_similarity_margin': ordered[:, -1] - ordered[:, -2],
                           'best_centroid_similarity': ordered[:, -1]})
    peer = z @ model['context_mean_ranks'].T
    for i, context in enumerate(CONTEXTS):
        result['rank_centroid_similarity_' + context] = similarities[:, i]
        result['mean_cell_spearman_' + context] = peer[:, i]
    for i in range(pcs.shape[1]):
        result[f'PC{i+1}'] = pcs[:, i]
    all_distances = np.empty((len(z), 3))
    for i, context in enumerate(CONTEXTS):
        delta = pcs - model['pc_context_means'][i]
        all_distances[:, i] = np.einsum('ij,ij->i', delta @ model['pc_context_cov_inverse'][i], delta)
        result['pc8_mahalanobis_squared_' + context] = all_distances[:, i]
    mahal = all_distances[np.arange(len(z)), labels]
    result['pc8_mahalanobis_squared'] = mahal
    result['outside_assigned_context_pc8_95pct'] = mahal > model['pc_context_mahal_p95'][labels]
    result['outside_all_context_pc8_95pct'] = (all_distances > model['pc_context_mahal_p95']).all(axis=1)
    result['below_assigned_context_similarity_p05'] = ordered[:, -1] < model['context_similarity_p05'][labels]
    return result


def make_plot(cells, external, pooled, explained, summary):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), layout='constrained')
    order = np.random.default_rng(20261004).permutation(len(cells))
    shuffled = cells.iloc[order]
    points = pd.concat([cells[['PC1','PC2']], external[['PC1','PC2']], pooled[['PC1','PC2']]])
    padding = (points.max() - points.min()) * .05
    for ax in axes.flat:
        for context, color in zip(CONTEXTS, COLORS):
            sub = shuffled[shuffled.context == context]
            ax.scatter(sub.PC1, sub.PC2, c=color, s=1, alpha=.12, rasterized=True, label=context)
        anchors = pooled[pooled.profile_kind == 'competition_pooled']
        for context, color in zip(CONTEXTS, COLORS):
            row = anchors[anchors.profile_name == 'context_' + context].iloc[0]
            ax.scatter(row.PC1, row.PC2, c=color, marker='P', s=90, edgecolor='black', linewidth=.5, zorder=5)
        ax.set_xlim(points.PC1.min()-padding.PC1, points.PC1.max()+padding.PC1)
        ax.set_ylim(points.PC2.min()-padding.PC2, points.PC2.max()+padding.PC2)
        ax.set_xlabel(f'PC1 ({explained[0]:.1%})'); ax.set_ylabel(f'PC2 ({explained[1]:.1%})')
    axes[0,0].set_title('Competition controls and pooled anchors\n55,200 cells; plus markers = pooled A/B/C')
    axes[0,0].legend(markerscale=6, title='Provided context')
    for ax, name in zip(list(axes.flat)[1:], REFERENCES):
        sub = external[external.reference == name]
        ax.scatter(sub.PC1, sub.PC2, c='black', s=5, alpha=.35, rasterized=True, label='External units')
        row = pooled[pooled.profile_name == name].iloc[0]
        ax.scatter(row.PC1, row.PC2, marker='*', s=240, c='#bb3377', edgecolor='black', linewidth=.7, zorder=6)
        if name == 'CD4':
            for condition, marker in [('Rest','s'), ('Stim8hr','^'), ('Stim48hr','D')]:
                part = pooled[pooled.profile_name == 'CD4_' + condition].iloc[0]
                ax.scatter(part.PC1, part.PC2, marker=marker, c='#bb3377', s=45, edgecolor='black', linewidth=.5, zorder=6)
        row = summary[(summary.reference == name) & (summary.stratum == 'all')].iloc[0]
        fractions = '/'.join(f'{row["fraction_"+c]:.1%}' for c in CONTEXTS)
        kind = 'single cells' if name in ['H1','HCT116','HEK293T'] else 'pseudobulk units'
        ax.set_title(f'{name}: {len(sub):,} {kind}\nNearest A/B/C: {fractions}; star = pooled profile', fontsize=10)
    fig.suptitle('External controls projected into the competition-only rank space\n7,447 jointly measured genes; zeros retained; coordinates frozen before external mapping', fontsize=14)
    fig.savefig(OUT / 'external_rank_projection.png', dpi=180)
    fig.savefig(OUT / 'external_rank_projection.pdf')
    plt.close(fig)
    profiles = pooled[pooled.profile_kind == 'external_pooled'].copy()
    fig, axes = plt.subplots(1,2,figsize=(13,5),layout='constrained')
    values = profiles[['pooled_profile_spearman_'+c for c in CONTEXTS]].to_numpy()
    im = axes[0].imshow(values, cmap='viridis', vmin=0, vmax=1, aspect='auto')
    axes[0].set_xticks(range(3), CONTEXTS); axes[0].set_yticks(range(len(profiles)), profiles.profile_name)
    for i in range(len(profiles)):
        for j in range(3):
            axes[0].text(j,i,f'{values[i,j]:.3f}',ha='center',va='center',color='white' if values[i,j]<.6 else 'black')
    fig.colorbar(im,ax=axes[0],label='Raw pooled-profile Spearman rho')
    axes[0].set_title('Pooled controls versus pooled A/B/C\nSame genes and raw ranks')
    main = summary[summary.stratum == 'all'].set_index('reference').loc[REFERENCES]
    base = np.zeros(len(main))
    for context,color in zip(CONTEXTS,COLORS):
        axes[1].barh(REFERENCES,main['fraction_'+context],left=base,color=color,label=context)
        base += main['fraction_'+context].to_numpy()
    axes[1].set_xlim(0,1); axes[1].set_xlabel('Fraction of external units'); axes[1].legend(title='Nearest context')
    axes[1].set_title('Unit-level assignment to frozen rank centroids\nFractions count units, not represented cells')
    fig.savefig(OUT/'external_mapping_summary.png',dpi=180)
    fig.savefig(OUT/'external_mapping_summary.pdf')
    plt.close(fig)


def write_results(summary, pooled, report):
    main = summary[summary.stratum == 'all'][['reference','units','unit_kind','fraction_A','fraction_B','fraction_C','median_similarity_margin','outside_all_pc8_fraction']]
    pool = pooled[pooled.profile_kind == 'external_pooled'][['profile_name','nearest_context','nearest_pooled_context'] + ['pooled_profile_spearman_'+c for c in CONTEXTS]]
    text = '# External controls in the A/B/C rank space\n\n'
    text += 'A joint, comparable space uses the **7,447 genes measured by all five references and A/B/C**. No missing genes were zero-filled. The original 18,533-gene coordinates cannot be projected exactly when external genes are unmeasured. All 55,200 competition cells were reranked on the common axis, then the cluster centroids and eight-component PCA were frozen before external mapping.\n\n'
    text += f"The competition-only shared-axis clustering has context ARI **{report['competition_context_ARI']:.6f}**, and independent-seed repeat ARI **{report['repeat_ARI']:.6f}**. The external references do not influence the fit.\n\n"
    text += '## Where external control units map\n\n' + table(main) + '\n\n'
    text += 'Fractions count external observations. HCT116, HEK293T and H1 observations are single cells. K562 observations are reconstructed pseudobulk groups; CD4 observations are raw pseudobulk sums. Their 75,094 and 34,620 represented cells are not available as separate cell rows in these releases. Assignments mean closest among A/B/C; they are not known cell identities or calibrated probabilities.\n\n'
    text += '## Pooled-reference comparison\n\n' + table(pool) + '\n\n'
    text += '`nearest_context` uses the frozen unit centroid on rank features. `nearest_pooled_context` independently compares raw pooled expression profiles with raw pooled A/B/C profiles using ordinary Spearman rho. Pooling and single-cell rank averaging are different operations, so they can disagree.\n\n'
    text += '## CD4 condition sensitivity\n\n' + table(summary[summary.stratum != 'all'][['reference','stratum','units','fraction_A','fraction_B','fraction_C','median_similarity_margin','outside_all_pc8_fraction']]) + '\n\n'
    text += '## Interpretation and distribution shift\n\nThe plot retains external positions even when they lie outside all competition clouds. Eight-PC Mahalanobis distances are compared with the empirical 95th percentile of the assigned competition context. The report also compares each external observation with all three PC clouds; `outside_all_pc8_fraction` flags those beyond every empirical 95th-percentile region. These fractions are diagnostics of displacement from supplied controls; it is not a biological identity test. Pseudobulks have much greater count depth than individual competition cells and can be displaced simply by pooling, while batch/platform/state differences add further shifts.\n\n'
    text += 'Row-level files include each external unit, its guide/batch/condition metadata, three rank-centroid similarities, three exact mean Spearman correlations to every competition context, margins, PCs and distribution-shift flags. All comparisons preserve measured zeros, use average tied ranks, and use no pseudocounts, logs, CPM or gene imputation.\n\n'
    text += '![External reference mapping](analysis/external_rank_mapping/external_rank_projection.png)\n\n![Mapping summary](analysis/external_rank_mapping/external_mapping_summary.png)\n\nSee [EXTERNAL_RANK_MAPPING.md](EXTERNAL_RANK_MAPPING.md) for reproduction, provenance and verification.\n'
    (ROOT / 'EXTERNAL_RANK_MAPPING_RESULTS.md').write_text(text, encoding='utf-8')


def run(seed=20261004, threads=12):
    OUT.mkdir(parents=True, exist_ok=True)
    axis, profiles, records = load_profiles()
    extraction = json.loads((CACHE/'extraction.json').read_text())
    for name, sha in extraction['outputs'].items():
        if digest(CACHE/name) != sha:
            raise ValueError('Extracted control cache changed')
    np.testing.assert_array_equal(pd.read_csv(CACHE/'gene_axis.csv').gene_name, axis)
    external_counts = sparse.load_npz(CACHE/'counts.npz')
    external_meta = pd.read_csv(CACHE/'units.csv')
    pd.DataFrame({'gene_name':axis}).to_csv(OUT/'gene_axis.csv',index=False)
    genes = pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    positions = {g:i for i,g in enumerate(genes)}
    axis_positions = np.array([positions[g] for g in axis])
    z = np.empty((55200,len(axis)),dtype=np.float32)
    cells = []; pooled_counts = []; inputs = []; offset=0
    original = json.loads((ROOT/'analysis/full_rank_clustering/analysis_report.json').read_text())
    for context in CONTEXTS:
        path = ROOT / f'context_{context}.h5ad'
        sha = digest(path)
        if sha != next(x['sha256'] for x in original['inputs'] if x['file'] == path.name):
            raise ValueError('Competition input checksum changed')
        obj = ad.read_h5ad(path,backed='r')
        try:
            if obj.shape != (18400,18533) or not np.array_equal(obj.var_names,genes) or set(obj.obs.context.astype(str)) != {context} or set(obj.obs.target_gene.astype(str)) != {'non-targeting'}:
                raise ValueError('Competition input annotations or axes differ')
            meta = obj.obs.copy(); meta.insert(0,'cell_id',obj.obs_names)
            totals=[];detected=[]; count_sum=np.zeros(len(axis))
            for start in range(0,len(obj),512):
                end=min(start+512,len(obj)); raw=sparse.csr_matrix(obj.X[start:end],dtype=np.float64)
                if not np.isfinite(raw.data).all() or (raw.data<0).any() or (raw.data!=np.floor(raw.data)).any():
                    raise ValueError('Invalid original counts')
                sub = raw[:,axis_positions].toarray()
                z[offset+start:offset+end]=rank_unit_rows(sub)
                count_sum += sub.sum(axis=0)
                totals.extend(sub.sum(axis=1));detected.extend((sub>0).sum(axis=1))
            meta['common_total_counts']=totals;meta['common_detected_genes']=detected
            cells.append(meta.reset_index(drop=True));pooled_counts.append(count_sum)
            inputs.append({'file':path.name,'sha256':sha});offset+=len(obj)
            print('Prepared competition rank rows:',context,offset,flush=True)
        finally:
            obj.file.close()
    cells=pd.concat(cells,ignore_index=True)
    with threadpool_limits(limits=threads):
        labels, centers, history, converged=cluster(z,3,seed)
        repeated, _, repeat_history, repeat_converged=cluster(z,3,seed+1)
        if not converged or not repeat_converged:
            raise ValueError('Shared-axis spherical fit did not converge')
        counts=pd.crosstab(cells.context,labels).reindex(index=CONTEXTS,columns=range(3),fill_value=0)
        context_rows, cluster_cols=linear_sum_assignment(-counts.to_numpy())
        cluster_order=cluster_cols[np.argsort(context_rows)]
        centroids=centers[cluster_order]
        cluster_to_context={int(k):CONTEXTS[i] for i,k in enumerate(cluster_order)}
        cells['blind_cluster']=labels
        cells['nearest_context']=[cluster_to_context[int(k)] for k in labels]
        pca=PCA(n_components=8,svd_solver='randomized',iterated_power=3,n_oversamples=10,random_state=seed)
        pca.fit(z);pcs=pca.transform(z)
        mean_ranks=np.array([z[cells.context==c].mean(axis=0,dtype=np.float64) for c in CONTEXTS])
        pc_means=[];cov_inv=[];mahal_p95=[];sim_p05=[]
        for i,context in enumerate(CONTEXTS):
            member=cells.context.to_numpy()==context
            pcmean=pcs[member].mean(axis=0,dtype=np.float64)
            inverse=np.linalg.pinv(np.cov(pcs[member],rowvar=False),hermitian=True)
            delta=pcs[member]-pcmean
            distance=np.einsum('ij,ij->i',delta@inverse,delta)
            pc_means.append(pcmean);cov_inv.append(inverse);mahal_p95.append(np.quantile(distance,.95))
            sim_p05.append(np.quantile(z[member]@centroids[i],.05))
        model={'centroids':centroids,'context_mean_ranks':mean_ranks,
               'pca_mean':pca.mean_,'pca_components':pca.components_,
               'pc_context_means':np.array(pc_means),'pc_context_cov_inverse':np.array(cov_inv),
               'pc_context_mahal_p95':np.array(mahal_p95),'context_similarity_p05':np.array(sim_p05)}
        for i in range(8):cells[f'PC{i+1}']=pcs[:,i]
        for i,c in enumerate(CONTEXTS):cells['rank_centroid_similarity_'+c]=z@centroids[i]
        # External controls transform only, using the frozen competition-only model.
        frames=[]
        for start in range(0,len(external_meta),256):
            end=min(start+256,len(external_meta))
            frames.append(frozen_projection(external_counts[start:end].toarray(),model))
        external=pd.concat([external_meta.reset_index(drop=True),pd.concat(frames,ignore_index=True)],axis=1)
        # external_meta contains no pre-existing mapping columns.
        pooled_raw=np.array(pooled_counts+[profiles[n].loc[axis].to_numpy() for n in REFERENCES+['CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']])
        pooled=frozen_projection(pooled_raw,model)
        pooled.insert(0,'profile_name',['context_'+c for c in CONTEXTS]+REFERENCES+['CD4_Rest','CD4_Stim8hr','CD4_Stim48hr'])
        pooled.insert(1,'profile_kind',['competition_pooled']*3+['external_pooled']*8)
        rho=np.clip(rank_unit_rows(pooled_raw)@rank_unit_rows(np.array(pooled_counts)).T, -1, 1)
        for i,c in enumerate(CONTEXTS):pooled['pooled_profile_spearman_'+c]=rho[:,i]
        pooled['nearest_pooled_context']=np.asarray(CONTEXTS)[rho.argmax(axis=1)]
    if external.columns.duplicated().any() or external.unit_id.duplicated().any():
        raise ValueError('Duplicate external columns or identifiers')
    summaries=[]
    for name in REFERENCES:
        sub=external[external.reference==name]
        groups=[('all',sub)]
        if name=='CD4':groups.extend((str(c),part) for c,part in sub.groupby('culture_condition',observed=True))
        for stratum,part in groups:
            row={'reference':name,'stratum':stratum,'units':len(part),'unit_kind':part.unit_kind.iloc[0],
                 'represented_cells':int(part.n_cells.sum()),
                 'median_similarity_margin':float(part.centroid_similarity_margin.median()),
                 'median_best_similarity':float(part.best_centroid_similarity.median()),
                 'outside_pc8_fraction':float(part.outside_assigned_context_pc8_95pct.mean()),
                 'outside_all_pc8_fraction':float(part.outside_all_context_pc8_95pct.mean()),
                 'low_similarity_fraction':float(part.below_assigned_context_similarity_p05.mean())}
            for c in CONTEXTS:
                row['fraction_'+c]=float((part.nearest_context==c).mean())
                row['median_mean_cell_spearman_'+c]=float(part['mean_cell_spearman_'+c].median())
            summaries.append(row)
    summary=pd.DataFrame(summaries)
    composition=counts.rename_axis(index='context',columns='blind_cluster').stack().reset_index(name='cells')
    outputs={'competition_cells_shared_space.csv':cells,'external_control_units_mapped.csv':external,
             'external_mapping_summary.csv':summary,'pooled_profiles_mapped.csv':pooled,
             'competition_cluster_composition.csv':composition,
             'pca_variance.csv':pd.DataFrame({'PC':np.arange(1,9),'explained_variance_ratio':pca.explained_variance_ratio_}),
             'pooled_raw_counts.csv':pd.DataFrame({'gene_name':axis,**{n:pooled_raw[i] for i,n in enumerate(pooled.profile_name)}})}
    for name,frame in outputs.items():frame.to_csv(OUT/name,index=False)
    np.savez_compressed(OUT/'frozen_shared_rank_model.npz',**model)
    (OUT/'external_control_extraction.json').write_text(json.dumps(extraction,indent=2)+'\n',encoding='utf-8')
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'seed':seed,'threads':threads,
            'competition_cells':len(cells),'genes':len(axis),'external_units':len(external),'competition_inputs':inputs,
            'competition_context_ARI':float(adjusted_rand_score(cells.context,labels)),
            'repeat_ARI':float(adjusted_rand_score(labels,repeated)),
            'blind_cluster_to_context':cluster_to_context,'fit_scope':'Competition only; external rows never influence PCA or clusters',
            'gene_policy':'Intersection of measured symbols across all five references and competition; no imputation',
            'count_transform':'Average tied raw ranks, centered and unit normalized; no added constant/log/CPM',
            'source_gene_coverage':{n:int(np.isin(genes,profiles[n].index).sum()) for n in REFERENCES},
            'original_full_rank_model_sha256':digest(ROOT/'analysis/full_rank_clustering/cluster_models.npz'),
            'refinement_history':history,'repeat_refinement_history':repeat_history,
            'versions':{name:version(name) for name in ['numpy','scipy','pandas','anndata','scikit-learn','matplotlib','h5py','pyarrow']},
            'scripts':{name:digest(ROOT/name) for name in ['extract_external_mapping_controls.py','map_external_rank_controls.py','test_external_rank_mapping.py','validate_external_rank_mapping.py']}}
    make_plot(cells,external,pooled,pca.explained_variance_ratio_,summary)
    report['outputs']={p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in OUT.iterdir() if p.is_file() and p.name not in ['analysis_report.json','validation.json']}
    (OUT/'analysis_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    write_results(summary,pooled,report)
    print(summary[summary.stratum=='all'].to_string(index=False),flush=True)
    print('Completed frozen reference mapping',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261004)
    parser.add_argument('--threads',type=int,default=12)
    args=parser.parse_args()
    run(args.seed,args.threads)
