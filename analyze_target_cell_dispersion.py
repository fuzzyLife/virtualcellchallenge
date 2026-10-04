"""Cell-level rank similarity and dispersion on the requested 300 VCC targets."""
from pathlib import Path
from datetime import datetime, timezone
import json
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_target_rank_correlations import ROOT, REFERENCES, digest, pairwise_rank_rows, table

OUT = ROOT / 'analysis/target_cell_dispersion'
CONTEXTS = ['A', 'B', 'C']


def normalized_ranks(raw):
    """Finite common-axis ranks; zero rows denote undefined constant vectors."""
    raw = np.asarray(raw, dtype=np.float64)
    if raw.ndim != 2 or not np.isfinite(raw).all():
        raise ValueError('Moment calculation requires an entirely finite common gene axis')
    ranks = rankdata(raw, method='average', axis=1)
    ranks -= ranks.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranks, axis=1)
    defined = norms > 0
    ranks[defined] /= norms[defined, None]
    return ranks, defined


def peer_moments(z, defined, groups, batch_size=1024):
    """Exact per-cell mean and population SD of rho to every valid group peer.

    For unit rank rows, rho(i,j)=z_i dot z_j. Sum and second moment
    follow from sum(z_j) and sum(z_j outer z_j). Exclude a cell's self pair.
    No individual pair matrix or sampling is needed for these moments.
    """
    groups = np.asarray(groups)
    result = {}
    own_square = np.einsum('ij,ij->i', z, z)
    for group in sorted(set(groups)):
        member = (groups == group) & defined
        peers = z[member]
        count = int(member.sum()) - member.astype(int)
        sums = peers.sum(axis=0)
        second = peers.T @ peers
        first = z @ sums - member * own_square
        squares = np.empty(len(z))
        for start in range(0, len(z), batch_size):
            end = min(start + batch_size, len(z))
            squares[start:end] = np.einsum('ij,ij->i', z[start:end] @ second, z[start:end])
        squares -= member * own_square**2
        valid = defined & (count > 0)
        mean = np.full(len(z), np.nan)
        variance = np.full(len(z), np.nan)
        mean[valid] = first[valid] / count[valid]
        variance[valid] = np.maximum(squares[valid] / count[valid] - mean[valid]**2, 0)
        variance[valid & (count == 1)] = 0  # One peer has exactly zero dispersion.
        result[group] = (mean, np.sqrt(variance), count)
    return result


def distribution(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {'defined': 0}
    q = np.quantile(x, [.05, .25, .5, .75, .95])
    return {'defined': len(x), 'mean': x.mean(), 'variance': x.var(ddof=1) if len(x)>1 else np.nan,
            'sd': x.std(ddof=1) if len(x)>1 else np.nan, 'p05': q[0], 'p25': q[1],
            'median': q[2], 'p75': q[3], 'p95': q[4], 'iqr': q[3]-q[1],
            'mad': np.median(np.abs(x-q[2]))}


def row_pair_spearman(left, right):
    a, av = normalized_ranks(left)
    b, bv = normalized_ranks(right)
    rho = np.einsum('ij,ij->i', a, b)
    rho[~(av & bv)] = np.nan
    return np.clip(rho, -1, 1)


def make_plot(cells, pairs):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), layout='constrained')
    colors = ['#4477aa', '#ee7733', '#228833']
    data, positions, labels = [], [], []
    for i, context in enumerate(CONTEXTS):
        for j, name in enumerate(REFERENCES):
            data.append(cells.loc[cells.context==context, 'common231_spearman_'+name].dropna())
            positions.append(i*6+j)
            labels.append(context+'\n'+name)
    boxes=axes[0,0].boxplot(data, positions=positions, tick_labels=labels, showfliers=False, patch_artist=True)
    for i, patch in enumerate(boxes['boxes']):
        patch.set_facecolor(colors[i//5]); patch.set_alpha(.65)
    axes[0,0].set_title('Individual cells vs external pooled profiles\nSame 231 genes', fontsize=11)
    axes[0,0].set_ylabel('Spearman rho'); axes[0,0].tick_params(axis='x', labelsize=8)
    plt.setp(axes[0,0].get_xticklabels(), rotation=40, ha='right')
    pair_names=list(pairs.pair_group.unique())
    axes[0,1].boxplot([pairs.loc[pairs.pair_group==g,'spearman'].dropna() for g in pair_names], tick_labels=pair_names, showfliers=False)
    axes[0,1].set_title('Individual cell pairs: sampled distribution\n300 genes', fontsize=11)
    axes[0,1].set_ylabel('Spearman rho')
    for context, color in zip(CONTEXTS, colors):
        delta=cells.loc[cells.context==context, 'delta_CD4_minus_H1_300']
        axes[1,0].hist(delta, bins=np.linspace(-.4,.4,101), density=True, histtype='step', label=context, color=color)
    axes[1,0].axvline(0, color='black', linestyle='--', linewidth=1)
    axes[1,0].set_title('Per-cell CD4 minus H1 correlation\n300 genes', fontsize=11)
    axes[1,0].set_xlabel('Positive values favor CD4'); axes[1,0].set_ylabel('Density'); axes[1,0].legend(title='Context')
    data=[]; labels=[]
    for context in CONTEXTS:
        for peer in CONTEXTS:
            data.append(cells.loc[cells.context==context,'mean_cell_correlation_to_'+peer].dropna())
            labels.append(context+' to '+peer)
    axes[1,1].boxplot(data, tick_labels=labels, showfliers=False)
    axes[1,1].set_title('Each cell mean correlation to all cells in a context\n300 genes', fontsize=11)
    axes[1,1].set_ylabel('Mean Spearman rho'); axes[1,1].tick_params(axis='x', labelsize=8)
    fig.suptitle('Cell heterogeneity on the VCC target-gene panel\nBoxes: median and IQR; whiskers: 1.5 IQR; raw ranks, zeros retained', fontsize=14)
    fig.savefig(OUT/'cell_dispersion.png', dpi=180)
    fig.savefig(OUT/'cell_dispersion.pdf')
    plt.close(fig)


def run(sample_pairs=20000, seed=20261004):
    OUT.mkdir(parents=True, exist_ok=True)
    old=ROOT/'analysis/target_rank_correlations'
    old_report=json.loads((old/'analysis_report.json').read_text())
    for name in ['per_cell_target_rank_correlations.csv','target_reference_raw_profiles.csv']:
        if digest(old/name)!=old_report['outputs'][name]['sha256']:
            raise ValueError('Previous verified correlation output changed: '+name)
    cells=pd.read_csv(old/'per_cell_target_rank_correlations.csv')
    refs=pd.read_csv(old/'target_reference_raw_profiles.csv')
    targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    if not np.array_equal(refs.gene_name,targets) or len(targets)!=300:
        raise ValueError('Target axes differ')
    common=np.isfinite(refs[REFERENCES].to_numpy()).all(axis=1)
    if common.sum()!=231:
        raise ValueError('Expected 231 shared target genes')
    pd.DataFrame({'gene_name':targets,'in_common_reference_axis':common}).to_csv(OUT/'gene_axis.csv',index=False)
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    pos={gene:i for i,gene in enumerate(genes)}
    axis=np.array([pos[gene] for gene in targets])
    arrays=[]; gene_stats=[]; inputs=[]
    for context in CONTEXTS:
        path=ROOT/f'context_{context}.h5ad'
        expected=next(x for x in old_report['competition_inputs'] if x['file']==path.name)
        if digest(path)!=expected['sha256']:
            raise ValueError('Original competition input changed')
        obj=ad.read_h5ad(path, backed='r')
        sub=cells[cells.context==context]
        if not np.array_equal(obj.var_names.to_numpy(),genes) or not np.array_equal(obj.obs_names.to_numpy(),sub.cell_id):
            raise ValueError('Original axes differ')
        raw=np.empty((obj.n_obs,300))
        for start in range(0,obj.n_obs,1024):
            end=min(start+1024,obj.n_obs)
            raw[start:end]=sparse.csr_matrix(obj.X[start:end],dtype=np.float64)[:,axis].toarray()
        if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any():
            raise ValueError('Expected finite nonnegative integer raw counts')
        arrays.append(raw)
        total=raw.sum(axis=0)
        mean=raw.mean(axis=0); variance=raw.var(axis=0,ddof=1)
        fano=np.full(300,np.nan); np.divide(variance,mean,out=fano,where=mean>0)
        gene_stats.append(pd.DataFrame({'context':context,'gene_name':targets,'cells':len(raw),
            'mean_raw_count':mean,'sample_variance_raw_count':variance,'sd_raw_count':np.sqrt(variance),
            'fano_raw_counts':fano,'zero_fraction':(raw==0).mean(axis=0)}))
        loo=np.empty(len(raw))
        for start in range(0,len(raw),1024):
            end=min(start+1024,len(raw))
            loo[start:end]=row_pair_spearman(raw[start:end],total[None,:]-raw[start:end])
        cells.loc[sub.index,'spearman_own_context_leave_one_out_300']=loo
        for name in REFERENCES:
            score,n,status=pairwise_rank_rows(raw[:,common],refs.loc[common,name].to_numpy())
            cells.loc[sub.index,'common231_spearman_'+name]=score
        obj.file.close(); inputs.append(expected)
        print(context,'loaded and leave-one-out correlations completed',flush=True)
    raw=np.vstack(arrays)
    z,defined=normalized_ranks(raw)
    moments=peer_moments(z,defined,cells.context.to_numpy())
    for context,(mean,sd,count) in moments.items():
        cells['mean_cell_correlation_to_'+context]=mean
        cells['sd_cell_correlation_to_'+context]=sd
        cells['valid_peer_cells_in_'+context]=np.where(defined,count,0)
    cells['delta_CD4_minus_H1_300']=cells.spearman_CD4-cells.spearman_H1
    cols=['common231_spearman_'+name for name in REFERENCES]
    scores=cells[cols].to_numpy()
    safe_scores=np.where(np.isfinite(scores),scores,-np.inf)
    best=safe_scores.max(axis=1)
    ties=((safe_scores==best[:,None]) & np.isfinite(scores)).sum(axis=1)
    cells['highest_common231_reference']=[REFERENCES[i] if t==1 else ('tie' if t>1 else 'undefined') for i,t in zip(safe_scores.argmax(axis=1),ties)]
    top=np.sort(scores,axis=1)
    cells['common231_top_minus_second']=top[:,-1]-top[:,-2]
    cells.to_csv(OUT/'per_cell_dispersion.csv',index=False)
    pd.concat(gene_stats,ignore_index=True).to_csv(OUT/'target_gene_count_dispersion.csv',index=False)
    summaries=[]; technical=[]; winners=[]
    for context,sub in cells.groupby('context',sort=True):
        for axis_name, prefix in [('pairwise_available_targets','spearman_'),('common_231_targets','common231_spearman_')]:
            for name in REFERENCES:
                summaries.append({'context':context,'comparison':name,'gene_axis':axis_name,**distribution(sub[prefix+name])})
        for peer in CONTEXTS:
            summaries.append({'context':context,'comparison':'mean_to_all_cells_'+peer,'gene_axis':'300_targets',**distribution(sub['mean_cell_correlation_to_'+peer])})
        for name in ['spearman_own_context_leave_one_out_300','delta_CD4_minus_H1_300']:
            summaries.append({'context':context,'comparison':name,'gene_axis':'300_targets',**distribution(sub[name])})
        for metric in ['spearman_H1','spearman_CD4','delta_CD4_minus_H1_300','spearman_own_context_leave_one_out_300']:
            for covariate in ['total_counts','nonzero_target_genes']:
                valid=np.isfinite(sub[metric]) & np.isfinite(sub[covariate])
                rho=spearmanr(sub.loc[valid,metric],sub.loc[valid,covariate]).statistic
                technical.append({'context':context,'correlation_metric':metric,'covariate':covariate,'cells':int(valid.sum()),'spearman':rho})
        for ref in REFERENCES+['tie','undefined']:
            n=int((sub.highest_common231_reference==ref).sum())
            winners.append({'context':context,'reference':ref,'cells':n,'fraction':n/len(sub)})
    summary=pd.DataFrame(summaries)
    summary.to_csv(OUT/'correlation_dispersion_summary.csv',index=False)
    pd.DataFrame(technical).to_csv(OUT/'count_depth_associations.csv',index=False)
    pd.DataFrame(winners).to_csv(OUT/'highest_reference_fractions.csv',index=False)
    guides=cells.groupby(['context','ntc_id']).agg(cells=('cell_id','size'),median_delta=('delta_CD4_minus_H1_300','median'),
        median_CD4=('spearman_CD4','median'),median_H1=('spearman_H1','median'),
        median_own_context=('spearman_own_context_leave_one_out_300','median'),median_depth=('total_counts','median'))
    guides.to_csv(OUT/'guide_dispersion.csv')
    rng=np.random.default_rng(seed)
    pair_frames=[]; pair_summaries=[]
    for i,left in enumerate(CONTEXTS):
        for right in CONTEXTS[i:]:
            li=np.flatnonzero((cells.context.to_numpy()==left)&defined)
            ri=np.flatnonzero((cells.context.to_numpy()==right)&defined)
            l=rng.choice(li,sample_pairs,replace=True); r=rng.choice(ri,sample_pairs,replace=True)
            if left==right:
                same=l==r
                while same.any():
                    r[same]=rng.choice(ri,int(same.sum()),replace=True); same=l==r
                l,r=np.minimum(l,r),np.maximum(l,r)
            rho=np.einsum('ij,ij->i',z[l],z[r])
            group=left+'-'+right
            pair_frames.append(pd.DataFrame({'pair_group':group,'cell_context_left':left,'cell_id_left':cells.cell_id.to_numpy()[l],
                'cell_context_right':right,'cell_id_right':cells.cell_id.to_numpy()[r],'spearman':rho,'valid_gene_pairs':300}))
            mean,sd,n=moments[right]
            exact_mean=np.average(mean[li],weights=n[li])
            exact_second=np.average(sd[li]**2+mean[li]**2,weights=n[li])
            pair_summaries.append({'pair_group':group,'all_unique_valid_cell_pairs':int(n[li].sum())//(2 if left==right else 1),
                'all_pairs_mean':exact_mean,'all_pairs_population_sd':np.sqrt(max(exact_second-exact_mean**2,0)),
                'sampled_pairs_with_replacement':sample_pairs,'sample_unique_pairs':len(set(zip(l.tolist(),r.tolist()))),
                **{'sample_'+k:v for k,v in distribution(rho).items()}})
    pairs=pd.concat(pair_frames,ignore_index=True)
    pairs.to_csv(OUT/'sampled_individual_cell_pair_correlations.csv',index=False)
    pair_summary=pd.DataFrame(pair_summaries)
    pair_summary.to_csv(OUT/'cell_pair_dispersion_summary.csv',index=False)
    make_plot(cells,pairs)
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'cells':len(cells),'target_genes':300,'common_reference_genes':int(common.sum()),
        'method':'Pairwise-valid Spearman on raw counts, average ties; measured zeros retained; no transforms',
        'cell_peer_moments':'Exact first/second moments over every valid peer cell, excluding self; 300 genes; population SD',
        'cell_pair_quantiles':'Uniform pairs with replacement, self excluded, same-context unordered; sampled quantiles only',
        'sample_pairs_per_group':sample_pairs,'seed':seed,'constant_cell_vectors':int((~defined).sum()),
        'external_scope':'Previously prepared pooled profiles, not individual external cells',
        'raw_count_variance':'Sample variance ddof=1; depth-sensitive raw-count Fano factor, not biological overdispersion estimate',
        'inputs':inputs,'target_list_sha256':digest(ROOT/'pert_counts.csv'),'previous_report_sha256':digest(old/'analysis_report.json'),
        'script_sha256':digest(Path(__file__)),'outputs':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(OUT.iterdir()) if p.suffix in ['.csv','.png','.pdf']}}
    (OUT/'analysis_report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    winners=pd.DataFrame(winners)
    text='# Target-gene cell dispersion results\n\n'+f"All {len(cells):,} cells; exactly the 300 competition target genes. No added constant or transformations.\n\n"
    text+='## Per-cell external-profile correlation spread\n\n'+table(summary[summary.gene_axis=='pairwise_available_targets'][['context','comparison','mean','sd','p05','median','p95','iqr']])+'\n\n'
    text+='## Highest correlation on the same 231 genes\n\n'+table(winners.pivot(index='context',columns='reference',values='fraction').reindex(columns=REFERENCES).reset_index())+'\n\n'
    text+='These fractions compare basal expression associations, not cell-type assignments.\n\n'
    text+='## Individual cell-pair similarity\n\n'+table(pair_summary[['pair_group','all_unique_valid_cell_pairs','all_pairs_mean','all_pairs_population_sd','sample_p05','sample_median','sample_p95']])+'\n\n'
    text+='Means and population SDs use ALL valid cell pairs exactly; p05/median/p95 use the reproducible uniform sample of '+str(sample_pairs)+' pairs per group. Same-cell comparisons are excluded.\n\n'
    text+='## Cell-to-own-context profile\n\n'+table(summary[summary.comparison=='spearman_own_context_leave_one_out_300'][['context','mean','sd','p05','median','p95']])+'\n\n'
    text+='Each pooled profile excludes the cell being compared. Correlation to a pool and correlation to another individual cell have different noise levels.\n\n'
    text+='## Count-depth associations\n\n'+table(pd.DataFrame(technical))+'\n\n'
    text+='These associations show sensitivity to count depth and detected targets; they do not separate technical and biological variation. Guide-level summaries are also supplied.\n\n'
    text+='![Dispersion figure](analysis/target_cell_dispersion/cell_dispersion.png)\n\nSee [TARGET_CELL_DISPERSION.md](TARGET_CELL_DISPERSION.md) for method, schemas and reproduction commands. Every cell has a row in `per_cell_dispersion.csv`; individual sampled pairs retain both cell IDs.\n'
    (ROOT/'TARGET_CELL_DISPERSION_RESULTS.md').write_text(text,encoding='utf-8')
    print(pair_summary.to_string(index=False),flush=True)
    print(winners.to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-pairs',type=int,default=20000)
    parser.add_argument('--seed',type=int,default=20261004)
    args=parser.parse_args()
    if args.sample_pairs<1:
        parser.error('sample-pairs must be positive')
    run(args.sample_pairs,args.seed)
