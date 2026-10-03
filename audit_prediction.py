"""Audit the first submission and paired count-generation ablations.

This is a diagnostic against provided controls, not a hidden-ground-truth score.
Original prediction files and analysis records are never overwritten.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
import build_submission as build

ROOT = Path(__file__).resolve().parent


def normalized_counts(matrix):
    matrix = matrix.astype(np.float64).tocsr()
    libraries = np.asarray(matrix.sum(axis=1)).ravel()
    if (libraries <= 0).any():
        raise ValueError('Cannot normalize an empty cell')
    return matrix.multiply((1e6/libraries)[:, None]).tocsr()


class ControlRankTest:
    """Cached two-sided asymptotic Mann-Whitney U, with ties and continuity.

    Compress control value frequencies once. Each prediction is compared with
    the complete control distribution, without sorting controls again. This
    gives the same p-values as scipy.stats.mannwhitneyu(method='asymptotic').
    """
    def __init__(self, controls):
        controls = sparse.csc_matrix(controls, dtype=np.float64)
        self.n = controls.shape[0]
        self.columns = []
        self.means = np.asarray(controls.mean(axis=0)).ravel()
        for column in range(controls.shape[1]):
            values = controls.data[controls.indptr[column]:controls.indptr[column+1]]
            values, counts = np.unique(values, return_counts=True)
            zeros = self.n-int(counts.sum())
            if zeros:
                if len(values) and values[0] == 0:
                    counts[0] += zeros
                else:
                    values = np.r_[0., values]
                    counts = np.r_[zeros, counts]
            counts = counts.astype(np.float64)
            prefix = np.r_[0., counts.cumsum()]
            ties = float((counts**3-counts).sum())
            self.columns.append((values, prefix, ties))

    def test(self, predicted):
        predicted = sparse.csc_matrix(predicted, dtype=np.float64)
        if predicted.shape[1] != len(self.columns):
            raise ValueError('Prediction/control axes differ')
        m = predicted.shape[0]
        total = m+self.n
        zvalues = np.full(len(self.columns), -np.inf)
        means = np.asarray(predicted.mean(axis=0)).ravel()
        for column, (ctrl, prefix, control_ties) in enumerate(self.columns):
            pv = predicted.data[predicted.indptr[column]:predicted.indptr[column+1]]
            pv, pc = np.unique(pv, return_counts=True)
            zeros = m-int(pc.sum())
            if zeros:
                if len(pv) and pv[0] == 0:
                    pc[0] += zeros
                else:
                    pv, pc = np.r_[0., pv], np.r_[zeros, pc]
            pc = pc.astype(np.float64)
            left = np.searchsorted(ctrl, pv, side='left')
            right = np.searchsorted(ctrl, pv, side='right')
            control_equal = prefix[right]-prefix[left]
            u = float((pc*(prefix[left]+prefix[right])/2).sum())
            ties = control_ties + float((pc**3-pc+3*control_equal*pc*(control_equal+pc)).sum())
            variance = m*self.n/12*(total+1-ties/(total*(total-1)))
            if variance > 0:
                z = (abs(u-m*self.n/2)-.5)/np.sqrt(variance)
                zvalues[column] = z
        pvalues = np.minimum(1., 2*stats.norm.sf(zvalues))
        adjusted = stats.false_discovery_control(pvalues, method='bh')
        lfc = np.log2((means+1e-9)/(self.means+1e-9))
        return pvalues, adjusted, lfc


def check_rank_kernel():
    rng = np.random.default_rng(97)
    controls = rng.integers(0, 7, (137, 19)).astype(float)
    predicted = rng.integers(0, 9, (43, 19)).astype(float)
    controls[:, 0] = predicted[:, 0] = 0
    controls[::2, 1] = .125
    predicted[::3, 1] = .125
    test = ControlRankTest(sparse.csr_matrix(controls))
    actual, _, _ = test.test(sparse.csr_matrix(predicted))
    expected = stats.mannwhitneyu(predicted, controls, axis=0,
                                 method='asymptotic', use_continuity=True).pvalue
    np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=1e-11)
    return {'passed': True, 'columns': 19, 'maximum_pvalue_difference':float(np.abs(actual-expected).max())}


def log(message):
    print(time.strftime('%H:%M:%S')+' '+message, flush=True)


def read_reference(root):
    with np.load(root/'analysis/reference_model_arrays.npz') as arr:
        return {name:arr[name].copy() for name in ('effects','available','basal','dest_positions','cell_lines')}


def select_targets(reference, targets, count):
    effects = reference['effects'][0]
    available = reference['available'][0]
    energy = np.sqrt(np.mean((.3*effects)**2, axis=1))
    covered = np.flatnonzero(available)
    ordering = covered[np.argsort(energy[covered])]
    ncovered = min(count-2, len(ordering))
    chosen = ordering[np.rint(np.linspace(0, len(ordering)-1, ncovered)).astype(int)].tolist()
    missing = np.flatnonzero(~available)
    chosen.extend(missing[:count-len(chosen)].tolist())
    return chosen


def record_api_status(output):
    result = subprocess.run([sys.executable,'-m','vcc','status','QQnFEBGlOk2pFFwkgtJL','--json'],
                            capture_output=True,text=True,encoding='utf-8',timeout=60)
    if result.returncode:
        return {'retrieved':False,'reason':'Read-only CLI status query failed; see user-provided original score record.'}
    status = json.loads(result.stdout)
    build.write_json(output/'published_api_status.json',status)
    return {'retrieved':True,'status':status}


def source_summary(root, ref, genes, targets):
    effect = ref['effects'][0]
    available = ref['available'][0]
    positions = ref['dest_positions']
    downstream = np.ones(effect.shape, dtype=bool)
    index = {gene:i for i,gene in enumerate(genes[positions])}
    for i,target in enumerate(targets):
        if target in index:
            downstream[i,index[target]] = False
    values = .3*effect[downstream]
    counts = ((np.abs(.3*effect)>.1)&downstream).sum(axis=1)
    eligible = {}
    for context in ['A','B','C']:
        frame = pd.read_csv(root/f'analysis/gene_statistics_{context}.csv')
        keep = frame['mean_cpm'].to_numpy()>5
        eligible[context] = {'tested_gene_count':int(keep.sum()),
                             'reference_covered_tested_genes':int(keep[positions].sum()),
                             'reference_covered_tested_fraction':float(keep[positions].sum()/keep.sum())}
    return {'targets_total':len(targets),'targets_covered':int(available.sum()),
            'gene_fraction_covered':len(positions)/len(genes),'gene_count_covered':len(positions),
            'source_context_correlations':json.loads((root/'analysis/context_transfer_diagnostics.json').read_text()),
            'identical_effect_multipliers_across_contexts':True,
            'downstream_pair_count':int(downstream.sum()),
            'nonzero_downstream_pair_fraction':float(np.count_nonzero(values)/len(values)),
            'downstream_pairs_above_absolute_log2':{str(v):int((np.abs(values)>v).sum()) for v in [.01,.05,.1,.2,.5]},
            'covered_targets_pairs_above_point1_quantiles':build.quantiles(counts[available]),
            'scored_axis_coverage':eligible}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--targets',type=int,default=6)
    parser.add_argument('--full-panel',action='store_true',help='Audit submitted DE counts for all 300 targets; omit ablations.')
    parser.add_argument('--skip-api',action='store_true')
    parser.add_argument('--contexts',nargs='+',choices=['A','B','C'],help='Audit selected contexts independently')
    args=parser.parse_args()
    if args.targets < 4:
        parser.error('--targets must be at least 4')
    root=args.root.resolve()
    output=root/'analysis/model_audit'
    output.mkdir(exist_ok=True)
    summary_name = 'full_panel_summary.json' if args.full_panel else 'audit_summary.json'
    csv_name = 'full_panel_de.csv' if args.full_panel else 'paired_ablations.csv'
    if args.contexts:
        suffix='_'+'_'.join(args.contexts)
        summary_name=summary_name.replace('.json',suffix+'.json')
        csv_name=csv_name.replace('.csv',suffix+'.csv')
    manifest,genes,targets=build.read_spec(root)
    ref=read_reference(root)
    indices=list(range(len(targets))) if args.full_panel else select_targets(ref,targets,args.targets)
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),
             'scope':'Provided controls and predicted cells only; no hidden challenge response labels.',
             'limitations':['Provided controls may differ from held-out scoring controls.',
                            'Predictions are compared with all provided controls, including sampled input cells; not an independent performance validation.',
                            'Increasing DE calls is not evidence of correct biology or an improved leaderboard score.',
                            'The small target selection is an effect-strength coverage diagnostic, not a random panel estimate.'],
             'rank_kernel_test':check_rank_kernel(),
             'source':source_summary(root,ref,genes,targets),
             'selected_targets':targets[indices].tolist(),
             'api_status':{'retrieved':False} if args.skip_api else record_api_status(output),
             'prediction_path':'prediction.h5ad',
             'input_package_sha256':json.loads((root/'analysis/leaderboard_result.json').read_text())['package_sha256']}
    model=(ref['effects'].astype(float),ref['available'],ref['basal'].astype(float),
           ref['dest_positions'],{'shrinkage':.3,'temperature':0.},ref['cell_lines'].tolist(),{})
    strong_model=(*model[:4],{'shrinkage':1.,'temperature':0.},*model[5:])
    rows=[]
    prediction=ad.read_h5ad(root/'prediction.h5ad',backed='r')
    for ci,context in enumerate(manifest['contexts']):
        if args.contexts and context not in args.contexts:
            continue
        log(f'Loading {context} controls and preparing complete control rank distributions')
        control=ad.read_h5ad(root/f'context_{context}.h5ad')
        assert np.array_equal(control.var_names.to_numpy(),genes)
        assert set(control.obs['context'].astype(str))=={context}
        x=control.X.astype(np.int32).tocsr()
        libraries=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel()
        guides=control.obs['ntc_id'].astype(str).to_numpy()
        cpm=normalized_counts(x)
        mean_cpm=np.asarray(cpm.mean(axis=0)).ravel()
        keep=np.flatnonzero(mean_cpm>5)
        ranks=ControlRankTest(cpm[:,keep])
        factors,_=build.context_factors(mean_cpm,model,genes,targets,.9)
        strong,_=build.context_factors(mean_cpm,strong_model,genes,targets,.9)
        panel_gene_mask=~np.isin(genes,targets)
        control_profile=np.log1p(5e4*np.asarray(x.sum(axis=0)).ravel()/x.sum())
        log(f'{context}: testing {len(keep)} genes across {len(indices)} targets')
        for ti in indices:
            rng=np.random.default_rng(np.random.SeedSequence([20261003,ci,ti]))
            sampled=build.guide_stratified_sample(guides,libraries,400,rng)
            offset=(ci*len(targets)+ti)*400
            observed=prediction[offset:offset+400,:]
            assert set(observed.obs['context'].astype(str))=={context}
            assert set(observed.obs['target_gene'].astype(str))=={targets[ti]}
            submitted=observed.X.tocsr()
            variants={'submitted':submitted}
            if not args.full_panel:
                target_only=np.ones(len(genes))
                target_only[np.flatnonzero(genes==targets[ti])[0]]=.1
                variants.update(sampled_controls=x[sampled],
                                target_only=build.simulate_counts(x,sampled,target_only,np.random.default_rng(20261003+ti)),
                                stronger_transfer=build.simulate_counts(x,sampled,strong[ti],np.random.default_rng(20261003+ti)),
                                zero_aware_original_transfer=build.simulate_counts_zero_aware(x,sampled,factors[ti],
                                    np.random.default_rng(20261003+ti),np.asarray(x.mean(axis=0)).ravel(),float(libraries.mean())))
            for variant,matrix in variants.items():
                normalized=normalized_counts(matrix)
                _,adjusted,lfc=ranks.test(normalized[:,keep])
                downstream=genes[keep]!=targets[ti]
                sig=(adjusted<.05)&downstream
                reference_axis=np.isin(keep,ref['dest_positions'])
                profile=np.log1p(5e4*np.asarray(matrix.sum(axis=0)).ravel()/matrix.sum())
                mask=genes!=targets[ti]
                rows.append({'context':context,'target_gene':targets[ti],'variant':variant,
                             'has_source_response':bool(ref['available'][0,ti]),
                             'tested_genes':int(downstream.sum()),'significant_downstream_genes':int(sig.sum()),
                             'significant_up':int((sig&(lfc>0)).sum()),'significant_down':int((sig&(lfc<0)).sum()),
                             'significant_outside_reference_axis':int((sig&~reference_axis).sum()),
                             'median_abs_observed_log2fc':float(np.median(np.abs(lfc[downstream]))),
                             'p95_abs_observed_log2fc':float(np.quantile(np.abs(lfc[downstream]),.95)),
                             'mean_library':float(matrix.sum()/matrix.shape[0]),
                             'log_profile_rms_to_controls_excluding_own_target':float(np.sqrt(np.mean((profile[mask]-control_profile[mask])**2))),
                             'log_profile_rms_to_controls_excluding_all_panel_targets':float(np.sqrt(np.mean((profile[panel_gene_mask]-control_profile[panel_gene_mask])**2)))})
            if not args.full_panel or (ti+1)%25==0:
                log(f'{context}/{targets[ti]}: DE calls '+str({r['variant']:r['significant_downstream_genes'] for r in rows if r['context']==context and r['target_gene']==targets[ti]}))
            pd.DataFrame(rows).to_csv(output/csv_name,index=False)
        del control,x,cpm,ranks,factors,strong
    prediction.file.close()
    frame=pd.DataFrame(rows)
    summary['group_summary']=frame.groupby(['context','variant'])['significant_downstream_genes'].agg(['count','min','median','mean','max']).reset_index().to_dict('records')
    summary['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    build.write_json(output/summary_name,summary)
    log(f'Audit complete: {output}')


if __name__=='__main__':
    main()
