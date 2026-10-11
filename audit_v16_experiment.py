"""Independent V16 identity, calibration, cache and exact null/baseline audit."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
from v16_methods import ROOT,A,FAMILIES,source_cache,predict
from run_v16_experiment import load_plan,digest,write
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v11_methods import emit as original_emit,SETTINGS
from evaluate_atlas_transfer import load_atlas_model

def main():
    plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];parent=ROOT/'experiments/v15_weighted_transport/analysis';old=json.loads((parent/'plan.json').read_text())
    model_audits=[];cache_audits=[];splits=[];baseline=0;null_identity=0;replayed=0;atlas,_=load_atlas_model()
    for f in FAMILIES:
        path=A/'models'/(str(f)+'.json')
        if not path.exists():continue
        r=json.loads(path.read_text());m=r['model'];assert r['inputs']['plan']==b.sha256(A/'plan.json') and m['excluded']==[str(f)] and m['reserved_label_targets_excluded']==plan['reserved_targets']
        n=np.sum([fold['numerator'] for fold in m['inner_folds']],axis=0);d=np.sum([fold['denominator'] for fold in m['inner_folds']],axis=0);ridge=.25*d.sum()/18
        expected=np.clip(np.divide(n+ridge,d+ridge,out=np.ones(18),where=d+ridge>0),0,1);expected[16:18]=1
        np.testing.assert_array_equal(n,m['numerator']);np.testing.assert_array_equal(d,m['denominator']);np.testing.assert_allclose(expected,m['gains'],rtol=1e-14,atol=1e-14)
        assert set(m['training_families'])==set(FAMILIES)-{str(f)}
        for fold in m['inner_folds']:assert set(fold['excluded_predictor_families'])=={str(f),fold['family']}
        model_audits.append({'held_family':str(f),'manual_gain_reconstruction':True,'outer_and_inner_exclusions_verified':True})
    for side in (A/'source_banks_cache').glob('*.json'):
        r=json.loads(side.read_text());file=side.with_suffix('.npz');assert r['sha256']==b.sha256(file) and r['inputs']['code']==b.sha256(ROOT/'v16_methods.py') and r['inputs']['source']==plan['source_cache_sha256'];assert set(r['training_families'])==set(FAMILIES)-set(r['inputs']['excluded'])
        with np.load(file,allow_pickle=False) as z:
            assert z['base'].shape==(300,18533) and np.isfinite(z['base']).all() and (z['votes']<=5).all() and ((z['agreement']>=0)&(z['agreement']<=1)).all();np.testing.assert_array_equal(z['base'][np.arange(300),z['own']],0.)
    for stage in ['null','development','confirmation']:
        report=A/f'{stage}_evaluation.json'
        if not report.exists():continue
        r=json.loads(report.read_text());assert r['complete'] and r['scores_sha256']==b.sha256(A/f'{stage}_scores.csv') and r['panels_sha256']==b.sha256(A/f'{stage}_panels.csv') and r['plan_sha256']==b.sha256(A/'plan.json')
        frame=pd.read_csv(A/f'{stage}_scores.csv');assert len(frame)==r['population_comparisons']
        for side in (A/(stage+'_cache')).glob('*.json'):
            rec=json.loads(side.read_text());file=side.with_suffix('.npz');assert rec['sha256']==b.sha256(file) and rec['inputs']['plan']==b.sha256(A/'plan.json')
            with np.load(file,allow_pickle=False) as z:assert all(z[k].shape==(18533,) and np.isfinite(z[k]).all() for k in z.files)
            cache_audits.append({'file':str(file.relative_to(ROOT)),'sha256':rec['sha256']})
    for ctx,prior in zip(plan['contexts'],old['contexts']):
        f=ctx['family'];assert f==prior['family'] and ctx['left']==prior['left'] and ctx['right']==prior['right'];x,meta,measured,_=raw_inputs(f,genes);assert meta.cell_id.is_unique and not set(ctx['left'])&set(ctx['right'])
        assert meta.iloc[ctx['left']+ctx['right']].gene_target.eq('Non-Targeting').all()
        for case in ctx['cells']:
            previous=next(c for c in prior['cells'] if c['target']==case['target'] and c['seed']==case['seed'] and c['split']==case['split']);assert all(previous[k]==case[k] for k in ['actual','template','reference'])
            assert len(set(case['actual']))==len(set(case['template']))==400 and not set(case['actual'])&set(ctx['left']+ctx['right']) and set(case['template'])<=set(ctx['left']);assert meta.iloc[case['actual']].gene_target.eq(case['target']).all()
            if case['split']=='development' and (A/'development_evaluation.json').exists():
                ti=int(np.flatnonzero(targets==case['target'])[0]);stem=f'{f}_{ti:03d}_{case["seed"]}';new=json.loads((A/'development_cache'/(stem+'.json')).read_text());old_rec=json.loads((parent/'development_cache'/(stem+'.json')).read_text())
                a=next(v for v in new['records'] if v['candidate']=='v13');c=next(v for v in old_rec['records'] if v['candidate']=='v13')
                for key in ['count_sha256','fidelity','reach','jaccard','nmae']:assert a[key]==c[key]
                with np.load(A/'development_cache'/(stem+'.npz'),allow_pickle=False) as z,np.load(parent/'development_cache'/(stem+'.npz'),allow_pickle=False) as q:
                    for key in ['truth','truth_jack','v13_pred','v13_jack']:np.testing.assert_array_equal(z[key],q[key])
                baseline+=1
        for case in ctx['null_cases']:
            assert len(set(case['template']))==400 and len(set(case['actual']))==200 and not set(case['template'])&set(case['actual']) and set(case['template']+case['actual'])<=set(ctx['left'])
            if (A/'null_evaluation.json').exists():
                rec=json.loads((A/'null_cache'/f'{f}_-01_{case["seed"]}.json').read_text());r=next(v for v in rec['records'] if v['candidate']=='balanced');template=sparse.csr_matrix(x[case['template']].toarray().astype(np.int64));assert digest(template)==r['count_sha256'];null_identity+=1
        if (A/'development_evaluation.json').exists():
            bank=source_cache(data,[f]);m=json.loads((A/'models'/(f+'.json')).read_text())['model'];_,p,bulk,_=moments(x[ctx['left']]);query=p*1e6;query[~measured]=np.nan;responses,_=predict(bank,m,query)
            first=ctx['development_targets'][0]
            for case in [v for v in ctx['cells'] if v['split']=='development' and v['target']==first]:
                ti=int(np.flatnonzero(targets==first)[0]);template=x[case['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);count,_=original_emit(template,np.arange(400),responses['v13'][ti],int(bank['own'][ti]),case['seed']+ti,SETTINGS['consensus20_noise65'],prepared=(libs,p,bulk),atlas=atlas)
                rec=json.loads((A/'development_cache'/f'{f}_{ti:03d}_{case["seed"]}.json').read_text());r=next(v for v in rec['records'] if v['candidate']=='v13');assert digest(count)==r['count_sha256'];replayed+=1
        splits.append({'family':f,'unique_source_cell_ids':True,'complete_control_halves_disjoint':True,'all_development_and_reserved_ids_match_parent_exactly':True,'all_treatments_templates_unique_and_disjoint':True});del x;gc.collect()
    if (A/'development_evaluation.json').exists():assert baseline==28 and replayed==6
    if (A/'null_evaluation.json').exists():assert null_identity==6
    prior_manifest=json.loads((parent/'report_manifest.json').read_text())
    for name,h in {**prior_manifest['inputs'],**prior_manifest['outputs']}.items():assert b.sha256(ROOT/name)==h,'Prior V15 evidence changed: '+name
    units=json.loads((A/'unit_tests.json').read_text());assert units['passed'] and units['tests_run']==9 and units['methods_sha256']==b.sha256(ROOT/'v16_methods.py')
    write(A/'independent_experiment_audit.json',{'passed':True,'models':model_audits,'splits':splits,'count_caches':cache_audits,'exact_reused_v15_baseline_count_profile_DE_cases':baseline,'original_baseline_same_input_regenerations':replayed,'bitwise_balanced_null_identity_cases':null_identity,'prior_v15_evidence_unchanged':True,'plan_sha256':b.sha256(A/'plan.json'),'script_sha256':b.sha256(Path(__file__)),'scope':'Original IDs, independent control splits, gain equations/exclusions, exact prior baselines/null identity, full cache hashes and selected baseline regeneration. Does not prove private accuracy.'})
    print('V16 independent audit passed:',baseline,'exact reused baselines;',replayed,'regenerated baselines;',null_identity,'exact null identities;',len(cache_audits),'caches',flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
