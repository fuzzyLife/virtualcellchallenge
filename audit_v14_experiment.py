"""Independent sample-coordinate ridge replay and complete V14 provenance checks."""
from pathlib import Path
import json,platform,sys
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v14_context_methods import ROOT,A,FAMILIES,family_profiles
from fit_v7_response import read_cache
from v7_response_methods import collapse
import build_submission as build


def main():
    plan=json.loads((A/'plan.json').read_text());execution=json.loads((A/'execution_freeze.json').read_text())
    for n,h in {**plan['dependencies'],**execution['dependencies']}.items():assert build.sha256(ROOT/n)==h,n
    history=json.loads((A/'metadata_correction.json').read_text());assert history['plan_sha256']==build.sha256(A/'plan.json')
    frozen=json.loads((A/'assessment_freeze.json').read_text());assert frozen['assessment_sha256']==build.sha256(ROOT/'assess_v14_experiment.py')
    old=json.loads((ROOT/'experiments/v13_response_programs/analysis/reference_model.json').read_text())
    for n,h in old['dependencies'].items():assert build.sha256(ROOT/n)==h,'V13 numerical dependency altered: '+n
    split_checks=[]
    for ctx in plan['contexts']:
        family=ctx['family'];folder=ROOT/('reference/v5_h1' if family=='H1' else f'reference/v7_sources/{family}')
        meta=pd.read_csv(folder/'selected_cells.csv')
        if family=='H1':
            meta['gene_target']=meta.target_gene.replace({'non-targeting':'Non-Targeting'})
            ids=meta.source_cell_id
        else:ids=meta.cell_id
        assert ids.is_unique,'Repeated original cell identities: '+family
        left=set(ctx['left']);right=set(ctx['right']);assert not left&right
        assert meta.iloc[sorted(left|right)].gene_target.eq('Non-Targeting').all()
        assert left|right==set(np.flatnonzero(meta.gene_target.eq('Non-Targeting')))
        assert not set(ctx['development_targets'])&set(plan['confirmation_targets'])
        for case in ctx['cells']:
            actual=case['actual'];template=case['template'];reference=case['reference']
            assert len(actual)==len(set(actual))==len(template)==len(set(template))==len(reference)==len(set(reference))==400
            assert meta.iloc[actual].gene_target.eq(case['target']).all()
            assert set(template)<=left and set(reference)<=right
            assert not set(actual)&(left|right)
        oracle=[c for c in plan['oracle_cases'] if c['family']==family]
        for case in oracle:
            calibration=set(case['calibration']);truth=set(case['truth'])
            assert len(calibration)==len(truth)==200 and not calibration&truth
            assert not (calibration|truth)&(left|right)
            assert meta.iloc[sorted(calibration|truth)].gene_target.eq(case['target']).all()
        split_checks.append({'family':family,'unique_source_cell_ids':len(ids),'control_halves_disjoint':True,
                             '400_cell_plan_cases_checked':len(ctx['cells']),'200_plus_200_oracle_cases_checked':len(oracle)})
    data=read_cache();profiles=family_profiles(data);own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
    adjusted=dict(data);adjusted['knockdown']=data['knockdown'].copy()
    source_own=data['profiles'][data['row_contexts'],own[data['row_targets']]]
    adjusted['knockdown'][~(np.isfinite(source_own)&(source_own>5))]=.9
    replays=[];inputs={}
    for family in ['H1','HCT116','HEK293T']:
        file=A/'models'/f'{family}.npz';side=file.with_suffix('.json');record=json.loads(side.read_text());assert build.sha256(file)==record['sha256']
        assert family not in record['metadata']['training_families']
        assert not (data['row_families'][record['metadata']['training_rows']]==family).any()
        with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
        train=np.flatnonzero(FAMILIES!=family);axis=arrays['f_axis'];assert not np.isin(data['genes'][axis],data['targets']).any()
        from scipy.stats import rankdata
        rank=np.array([rankdata(p[axis],method='average')/len(axis) for p in profiles[train]])
        np.testing.assert_allclose(arrays['f_center'],rank.mean(0),atol=1e-14)
        scores=((rank-arrays['f_center'])@arrays['f_basis'].T)/arrays['f_scale']
        np.testing.assert_allclose(np.mean(scores**2,axis=0),1.,atol=1e-10)
        y,var,valid,kd=collapse(adjusted,data['row_families']!=family)
        gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
        rng=np.random.default_rng(2026101431);targets=rng.choice(len(own),12,replace=False);columns=rng.choice(len(data['genes']),24,replace=False)
        for name,penalty,shuffle in [('context2',2.,False),('context8',8.,False),('shuffled2',2.,True)]:
            design=scores.copy()
            if shuffle:design=design[np.random.default_rng(2026101403).permutation(len(train))]
            max_error=0.
            for ti in targets:
                for gi in columns:
                    usable=valid[train,ti,gi];w=np.where(usable,np.clip(.04/np.maximum(var[train,ti,gi],.0025),.1,4.),0.)
                    if w.sum()==0:center=np.zeros(2);beta=np.zeros(2)
                    else:
                        w=w/w.sum();center=w@design;values=y[train,ti,gi];mean=w@values
                        xc=design-center;weights=w*usable.sum()
                        gram=xc.T@(weights[:,None]*xc)+penalty*np.eye(2)
                        beta=np.linalg.solve(gram,xc.T@(weights*(values-mean)))
                        if usable.sum()<3:beta[:]=0
                    if gi==own[ti]:beta[:]=0
                    error=float(np.max(abs(beta-arrays['s_'+name][ti,gi])));max_error=max(max_error,error)
                    np.testing.assert_allclose(arrays['c_'+name][ti,gi],center,atol=2e-6)
            assert max_error<2e-6,(family,name,max_error)
            np.testing.assert_array_equal(arrays['s_'+name][np.arange(len(own)),own],0.)
            replays.append({'family':family,'candidate':name,'sampled_target_gene_pairs':len(targets)*len(columns),'maximum_beta_difference':max_error})
        inputs[str(file.relative_to(ROOT))]=build.sha256(file)
        del arrays,y,var,valid
    caches=[]
    for directory in sorted(A.glob('*_cache')):
        for side in directory.glob('*.json'):
            d=json.loads(side.read_text());file=side.with_suffix('.npz');assert build.sha256(file)==d['sha256']
            with np.load(file,allow_pickle=False) as z:
                assert all(z[k].shape==(18533,) and np.isfinite(z[k]).all() for k in z.files)
            caches.append({'file':str(file.relative_to(ROOT)),'sha256':d['sha256']})
    aggregate_checks=[]
    from run_v14_experiment import FIELDS
    for stage,expected_cases,expected_rows,expected_panels in [('development',72,144,48),('oracle',27,27,9)]:
        evaluation=json.loads((A/f'{stage}_evaluation.json').read_text())
        assert evaluation['completed'] and evaluation['generated_or_real_population_cases']==expected_cases
        assert evaluation['scores_sha256']==build.sha256(A/f'{stage}_scores.csv')
        assert evaluation['summary_sha256']==build.sha256(A/f'{stage}_summary.csv')
        rows=pd.read_csv(A/f'{stage}_scores.csv');panels=pd.read_csv(A/f'{stage}_panels.csv')
        assert len(rows)==expected_rows and len(panels)==expected_panels
        primary=panels[panels.reference_mode=='full_independent']
        recomputed=primary.groupby(['candidate','family'])[FIELDS].mean().groupby('candidate').mean()
        saved=pd.read_csv(A/f'{stage}_summary.csv').set_index('candidate')[FIELDS]
        np.testing.assert_allclose(recomputed.loc[saved.index],saved,rtol=1e-12,atol=1e-12,equal_nan=True)
        aggregate_checks.append({'stage':stage,'population_cases':expected_cases,'score_rows':len(rows),'panel_rows':len(panels),'family_balanced_aggregation_replayed':True})
    selection=json.loads((A/'selection.json').read_text())
    assert selection['development_sha256']==build.sha256(A/'development_summary.csv')
    assert selection['panels_sha256']==build.sha256(A/'development_panels.csv')
    assert selection['promote']==all(selection['checks'].values())
    assert selection['metadata_correction_sha256']==build.sha256(A/'metadata_correction.json')
    if not selection['promote']:assert not (A/'confirmation_evaluation.json').exists()
    build.write_json(A/'independent_audit.json',{'passed':True,'ridge_coordinate_replays':replays,'cache_files_checked':len(caches),'split_identity_checks':split_checks,'aggregate_checks':aggregate_checks,
        'cache_hashes':caches,'model_hashes':inputs,'frozen_v13_dependencies_unchanged':True,'own_target_slopes_zero':True,
        'training_features_exclude_panel_genes_and_held_family':True,'metadata_correction_applied':True,
        'script_sha256':build.sha256(Path(__file__)),'runtime':{'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__},
        'scope':'All cached profile/jackknife arrays and SHA checks; independent scalar ridge replay at 2592 sampled coordinates. No claim to rerun every DE test or validate hidden biological accuracy.'})
    print('V14 independent audit passed:',len(caches),'count/profile caches;',len(replays)*288,'ridge coordinate replays')


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
