"""Freeze a bounded V6 model using target-held-out source generator validation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_context_ensemble import target_folds
from calibrate_v6_submission import candidates_summary
from v6_methods import CANDIDATES
ROOT=build.ROOT; A=ROOT/'experiments/v6_fidelity_calibrated/analysis'


def select_bounded(summary, permitted=None):
    baseline=summary[summary.candidate=='v4'].iloc[0]
    allowed=[k for k,v in CANDIDATES.items() if v['blend']<=.25 and (permitted is None or k in permitted)]
    eligible=summary[summary.candidate.isin(allowed)]
    eligible=eligible[(eligible.fidelity>=baseline.fidelity-.005)&(eligible.jaccard>=baseline.jaccard-.005)&
        (eligible.expression_ratio<=baseline.expression_ratio+1e-12)]
    if eligible.empty:return 'v4'
    score=eligible.fidelity+eligible.jaccard+.5*(1-eligible.expression_ratio)
    return str(eligible.iloc[int(np.argmax(score.to_numpy()))].candidate)


def assemble(old,new,setting):
    old=np.asarray(old,dtype=np.float64);new=np.asarray(new,dtype=np.float64)
    effect=old.copy()
    weight=setting['blend']
    if not 0<=weight<=.25:raise ValueError('V6 response blend exceeds prespecified bound')
    for ci in [1,2]:
        same=(old[ci]*new[ci]>0)
        effect[ci,same]=(1-weight)*old[ci,same]+weight*new[ci,same]
    return effect.astype(np.float32)


def main():
    scores=pd.read_csv(A/'source_generator_scores.csv').to_dict('records')
    profiles={}
    with np.load(A/'source_profile_cache.npz',allow_pickle=False) as z:
        genes,targets=z['genes'].copy(),z['targets'].copy()
        for candidate in CANDIDATES:
            for d in [0,1]:
                for s in [0,1]:
                    arrays=[z[f'{candidate}_{d}_{s}_{k}'] for k in ['pred','truth','control','pred_jack','truth_jack','control_jack']]
                    for i,target in enumerate(targets):profiles[(str(target),d,s,candidate)]=tuple(a[i] for a in arrays)
    _,axis,official_targets=build.read_spec(ROOT)
    # Reduced-size groups are reported, but settings are selected on groups with at least 200 distinct cells in both directions.
    size=pd.DataFrame(scores).groupby('target_gene').cells.min()
    primary=np.array([t for t in targets if size.loc[t]>=200])
    assert len(primary)>=5
    activity=pd.read_csv(A/'context_generator_ablation.csv')
    budget=activity.groupby(['context','candidate']).predicted_DE_calls.mean().unstack()
    generator_guards={k:bool((budget[k]/budget['v4']>=.9).all()) for k in ['v4_dual','v4_cal10','v4_cal25']}
    permitted=[k for k,v in CANDIDATES.items() if v['blend']<=.25 and generator_guards.get(k,True)]
    folds=target_folds(targets);outer=[];choices=[]
    for fold in range(5):
        all_train=targets[folds!=fold];all_test=targets[folds==fold]
        train=np.array([t for t in all_train if t in set(primary)]);test=np.array([t for t in all_test if t in set(primary)])
        assert len(train) and len(test)
        training=candidates_summary([r for r in scores if r['target_gene'] in set(train)],profiles,genes,train,official_targets)
        chosen=select_bounded(training,permitted)
        heldout=candidates_summary([r for r in scores if r['target_gene'] in set(test)],profiles,genes,test,official_targets)
        for row in heldout.to_dict('records'):outer.append({'fold':fold,'primary_targets':len(test),'selected':row['candidate']==chosen,**row})
        choices.append({'fold':fold,'candidate':chosen,'train_targets':all_train.tolist(),'test_targets':all_test.tolist(),'primary_train_targets':train.tolist(),'primary_test_targets':test.tolist()})
    outer=pd.DataFrame(outer);outer.to_csv(A/'bounded_outer_scores.csv',index=False)
    def weighted(frame):
        result={}
        for key in ['fidelity','reach','jaccard','nmae','expression_ratio','pds','uncorrected_profile_mse']:
            keep=frame[key].notna();result[key]=float(np.average(frame.loc[keep,key],weights=frame.loc[keep,'primary_targets'])) if keep.any() else None
        return result
    selected=weighted(outer[outer.selected]);baseline=weighted(outer[outer.candidate=='v4'])
    supported=selected['fidelity']>=baseline['fidelity']-.005 and selected['jaccard']>=baseline['jaccard']-.005 and selected['expression_ratio']<=baseline['expression_ratio']
    full=candidates_summary([r for r in scores if r['target_gene'] in set(primary)],profiles,genes,primary,official_targets)
    full.to_csv(A/'primary_source_candidate_summary.csv',index=False)
    chosen=select_bounded(full,permitted) if supported else 'v4';setting=CANDIDATES[chosen]
    neural=json.loads((A/'neural_template_selection.json').read_text())
    assert neural['production_neural_weight']==0
    oldpath=ROOT/'experiments/v4_lineage_response/analysis/reference_model_arrays.npz';newpath=ROOT/'experiments/v5_atlas_calibrated/analysis/reference_model_arrays.npz'
    oldmeta=json.loads(oldpath.with_name('reference_model.json').read_text());newmeta=json.loads(newpath.with_name('reference_model.json').read_text())
    assert build.sha256(oldpath)==oldmeta['model_arrays_sha256'] and build.sha256(newpath)==newmeta['model_arrays_sha256']
    with np.load(oldpath,allow_pickle=False) as z:
        old=z['effects'].copy();np.testing.assert_array_equal(z['genes'],axis);np.testing.assert_array_equal(z['targets'],official_targets)
    with np.load(newpath,allow_pickle=False) as z:
        new=z['effects'].copy();np.testing.assert_array_equal(z['genes'],axis);np.testing.assert_array_equal(z['targets'],official_targets)
    effect=assemble(old,new,setting)
    np.testing.assert_array_equal(effect[0],old[0])
    assert not ((effect*old<0)&(old!=0)).any() and (effect[old==0]==0).all()
    np.savez_compressed(A/'reference_model_arrays.npz',effects=effect,genes=axis.astype(str),targets=official_targets.astype(str))
    changes=[]
    for ci,context in enumerate('ABC'):
        for ti,target in enumerate(official_targets):
            delta=effect[ci,ti].astype(float)-old[ci,ti].astype(float)
            changes.append({'context':context,'target_gene':target,'changed_effects':int(np.count_nonzero(delta)),
                'maximum_absolute_change':float(np.abs(delta).max()),'opposite_signs_from_v4':int((effect[ci,ti]*old[ci,ti]<0).sum())})
    pd.DataFrame(changes).to_csv(A/'v4_v6_effect_changes.csv',index=False)
    result={'model':'V6 fidelity-guarded V4 response with source-validated count calibration',
        'candidate':chosen,'setting':setting,'context_generator_guards':generator_guards,'minimum_activity_fraction_of_v4':.9,'context_ablation_sha256':build.sha256(A/'context_generator_ablation.csv'),'nested_adoption_supported':bool(supported),'nested_procedure':selected,'nested_v4_comparator':baseline,
        'target_outer_choices':choices,'primary_source_targets':primary.tolist(),'primary_minimum_cells':200,'small_group_scope':'All smaller groups reported but not used to select settings.','neural_weight':0.,'neural_gate_sha256':build.sha256(A/'neural_template_selection.json'),
        'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'v4_model_sha256':build.sha256(oldpath.with_name('reference_model.json')),
        'v4_arrays_sha256':build.sha256(oldpath),'v5_arrays_sha256':build.sha256(newpath),'context_original_inputs':oldmeta['context_original_inputs'],
        'source_scores_sha256':build.sha256(A/'source_generator_scores.csv'),'source_cache_sha256':build.sha256(A/'source_profile_cache.npz'),
        'source_calibration_sha256':build.sha256(A/'source_selection.json'),'script_sha256':build.sha256(Path(__file__)),
        'methods_sha256':build.sha256(ROOT/'v6_methods.py'),'response_policy':'A unchanged; B/C blend only existing nonzero effects agreeing in sign, at most 25%. No promoter addition or newly nonzero V4 coordinates.',
        'limitations':['Source scores are unanchored H1 comparisons on 25 targets with reduced cell groups, not hidden A/B/C accuracy.',
            'The source H1 data and earlier validation results were previously inspected; this is a comparative ablation, not a new independent dataset.',
            'H1 donor/batch validation may not transfer to anonymous contexts; CD4 response policy is retained.',
            'Neural reconstruction is observational and was not converted into a causal intervention model.']}
    build.write_json(A/'reference_model.json',result)
    print('Frozen V6:',chosen,setting,'nested supported',supported,flush=True)
    print('Nested:',selected,'\nV4 baseline:',baseline,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
