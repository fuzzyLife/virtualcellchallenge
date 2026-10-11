"""Independent V15 provenance, split, source-quantile and coefficient-selection audit."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from v15_ensemble_methods import ROOT,A,FAMILIES,GRID
from fit_v7_response import read_cache,raw_inputs
from run_v15_experiment import load_plan

def sorted_quantile(a,q):
    a=np.sort(a,axis=0);v=q*(len(a)-1);lo=np.floor(v).astype(int);hi=np.ceil(v).astype(int);fraction=v-lo
    return (a[lo]*(1-fraction[:,None])+a[hi]*fraction[:,None]).T

def main():
    plan=load_plan();data=read_cache();own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);splits=[];quantiles=[];models=[];caches=[]
    for context in plan['contexts']:
        family=context['family'];x,meta,measured,provenance=raw_inputs(family,data['genes']);assert meta.cell_id.is_unique
        left=set(context['left']);right=set(context['right']);assert not left&right
        assert meta.iloc[sorted(left|right)].gene_target.eq('Non-Targeting').all()
        for case in context['cells']:
            for name in ['actual','template','reference']:assert len(case[name])==len(set(case[name]))==400
            assert set(case['template'])<=left and set(case['reference'])<=right and not set(case['actual'])&(left|right)
            assert meta.iloc[case['actual']].gene_target.eq(case['target']).all()
        for case in context['null_cases']:
            assert len(case['template'])==len(set(case['template']))==len(case['actual'])==len(set(case['actual']))==200
            assert set(case['template'])<=left and set(case['actual'])<=left and not set(case['template'])&set(case['actual'])
        assert not set(context['development_targets'])&set(plan['reserved_targets'])
        splits.append({'family':family,'source_cell_ids_unique':True,'control_halves_disjoint':True,'400_cell_cases':len(context['cells']),'null_200_plus_200_cases':len(context['null_cases'])})
        path=A/'distribution_sources'/(family+'.npz');side=json.loads(path.with_suffix('.json').read_text());assert side['sha256']==build.sha256(path)
        assert side['source']==provenance and side['inputs']['code']==build.sha256(ROOT/'v15_population_methods.py')
        with np.load(path,allow_pickle=False) as z:
            selected=np.flatnonzero(z['valid'].any(1))[:3];columns=np.random.default_rng(2026101543).choice(18533,24,replace=False)
            controls=np.flatnonzero(meta.gene_target.eq('Non-Targeting'));libs=np.asarray(x[controls].sum(1)).ravel();cp=np.asarray(x[controls][:,columns].multiply((1/libs)[:,None]).mean(0)).ravel();split=float(np.median(libs));assert split==float(z['depth_split'])
            qctrl=[];fi=int(np.flatnonzero(FAMILIES==family)[0]);maximum=0.
            for bi in range(2):
                available=controls[(libs<=split) if bi==0 else (libs>split)];ids=np.random.default_rng(2026101501+fi*10+bi).choice(available,min(len(available),800),replace=False)
                depths=np.asarray(x[ids].sum(1)).ravel();r=x[ids][:,columns].toarray()/depths[:,None];r=np.divide(r,cp[None,:],out=np.zeros_like(r),where=cp[None,:]>0)
                qctrl.append(sorted_quantile(r,np.array([0.,.25,.5,.75,1.])))
            for ti in selected:
                ids=np.flatnonzero(meta.gene_target.eq(data['targets'][ti]));depths=np.asarray(x[ids].sum(1)).ravel();r=x[ids][:,columns].toarray()/depths[:,None];p=r.mean(0);r=np.divide(r,p[None,:],out=np.zeros_like(r),where=p[None,:]>0)
                for bi in range(2):
                    use=(depths<=split) if bi==0 else (depths>split)
                    if use.sum()<25:continue
                    wanted=np.clip(sorted_quantile(r[use],np.array([0.,.25,.5,.75,1.]))-qctrl[bi],-4,4);wanted[~z['valid'][ti,columns]]=0
                    error=float(np.max(abs(wanted-z['delta'][ti,bi,columns])));maximum=max(maximum,error)
            assert maximum<2e-5,(family,maximum)
            quantiles.append({'family':family,'targets':len(selected),'genes':len(columns),'bins':2,'quantile_coordinates':len(selected)*len(columns)*2*5,'maximum_difference':maximum})
        del x;gc.collect()
    for side in (A/'models').glob('*.json'):
        d=json.loads(side.read_text());m=d['model'];assert d['inputs']['plan']==build.sha256(A/'plan.json');assert m['weights']==list(GRID[int(np.argmax(m['inner_mean_objectives']))])
        assert set(m['inner_families'])==set(FAMILIES)-set(m['excluded']);assert m['weights'][0]>=.5 and abs(sum(m['weights'])-1)<1e-12
        for excluded in m['excluded']:assert excluded not in m['inner_families']
        models.append({'family':side.stem,'weights':m['weights'],'excluded':m['excluded'],'grid_selection_replayed':True})
    for side in (A/'response_cache').glob('*.json'):
        rec=json.loads(side.read_text());file=side.with_suffix('.npz');assert build.sha256(file)==rec['sha256'];assert set(rec['training_families'])==set(FAMILIES)-set(rec['inputs']['excluded'])
        with np.load(file,allow_pickle=False) as z:
            for name in ['v13','student']:
                assert z[name].shape==(300,18533) and np.isfinite(z[name]).all();np.testing.assert_array_equal(z[name][np.arange(300),own],0.)
    for stage in ['null','development','confirmation']:
        report=A/f'{stage}_evaluation.json'
        if not report.exists():continue
        evaluation=json.loads(report.read_text());assert evaluation['complete']
        assert evaluation['scores_sha256']==build.sha256(A/f'{stage}_scores.csv') and evaluation['panels_sha256']==build.sha256(A/f'{stage}_panels.csv')
        scores=pd.read_csv(A/f'{stage}_scores.csv');assert len(scores)==evaluation['population_cases']
        for side in (A/(stage+'_cache')).glob('*.json'):
            rec=json.loads(side.read_text());file=side.with_suffix('.npz');assert rec['sha256']==build.sha256(file)
            with np.load(file,allow_pickle=False) as z:assert all(z[k].shape==(18533,) and np.isfinite(z[k]).all() for k in z.files)
            caches.append({'file':str(file.relative_to(ROOT)),'sha256':rec['sha256']})
    baseline_replays=0
    if (A/'development_evaluation.json').exists():
        old=pd.read_csv(ROOT/'experiments/v14_context_interactions/analysis/development_scores.csv')
        old=old[(old.candidate=='v13')&(old.reference_mode=='full_independent')]
        new=pd.read_csv(A/'development_scores.csv');new=new[new.candidate=='v13']
        pairs=old.merge(new,on=['family','target','seed'],suffixes=('_old','_new'))
        assert len(pairs)==18,'Expected all 18 shared V14/V15 baseline cases'
        assert (pairs.count_sha256_old==pairs.count_sha256_new).all()
        for field in ['fidelity','reach','jaccard','nmae']:
            np.testing.assert_allclose(pairs[field+'_old'],pairs[field+'_new'],atol=0,rtol=0,equal_nan=True)
        baseline_replays=len(pairs)
    build.write_json(A/'independent_experiment_audit.json',{'passed':True,'exact_v14_baseline_count_and_DE_replays':baseline_replays,'splits':splits,'quantile_replays':quantiles,'models':models,'count_profile_caches':caches,'v13_dependencies_unchanged':True,'script_sha256':build.sha256(Path(__file__)),
        'scope':'Original unique identities, frozen dependencies, all cache hashes/shapes, nested grid choices, manually sorted/interpolated source quantile samples. Does not prove hidden leaderboard accuracy or independently rerun every DE test.'})
    print('V15 experiment audit passed;',sum(v['quantile_coordinates'] for v in quantiles),'source quantile coordinates;',len(caches),'count caches',flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
