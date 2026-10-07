"""Validate V6 generators on disjoint public H1 batches and target-held-out settings."""
from pathlib import Path
import json
import gc
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from concurrent.futures import ThreadPoolExecutor
import build_submission as build
from evaluate_atlas_transfer import load_atlas_model
from fit_context_ensemble import target_folds
from v6_methods import CANDIDATES,moments,emit,de,direction_metrics,profile_jackknife,summarize_profiles

ROOT=build.ROOT; A=ROOT/'experiments/v6_fidelity_calibrated/analysis'; DATA=ROOT/'reference/v5_h1'
SEED=20261007


def candidates_summary(records, profiles, genes, targets, official_targets):
    frame=pd.DataFrame(records)
    rows=[]
    for candidate in CANDIDATES:
        group=frame[frame.candidate==candidate]
        metrics=group[['fidelity','reach','jaccard','nmae','predicted_calls','reference_calls']].mean().to_dict()
        panels=[]
        for direction in [0,1]:
            for seed in [0,1]:
                selected=[profiles[(str(t),direction,seed,candidate)] for t in targets]
                arr=[np.stack([r[k] for r in selected]) for k in range(6)]
                panels.append(summarize_profiles(*arr,genes,targets,official_targets))
        rows.append({'candidate':candidate,**metrics,**{k:float(np.mean([p[k] for p in panels])) for k in panels[0]}})
    return pd.DataFrame(rows)


def choose(summary):
    base=summary[summary.candidate=='v4'].iloc[0]
    # Fidelity is the primary adoption guard; source expression error must also improve.
    valid=summary[(summary.fidelity>=base.fidelity-.005)&(summary.jaccard>=base.jaccard-.005)&
        (summary.expression_ratio<=base.expression_ratio+1e-12)]
    score=valid.fidelity+valid.jaccard+.5*(1-valid.expression_ratio)
    order=np.argsort(-score.to_numpy(),kind='stable')
    return str(valid.iloc[order[0]].candidate)


def main():
    A.mkdir(parents=True,exist_ok=True)
    source=json.loads((DATA/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:
        assert build.sha256(DATA/name)==source[key]
    x=sparse.load_npz(DATA/'selected_counts.npz').astype(np.int32)
    meta=pd.read_csv(DATA/'selected_cells.csv'); genes=pd.read_csv(DATA/'gene_names.csv').gene_name.to_numpy(dtype=str)
    assert meta.source_cell_id.is_unique
    _,_,official_targets=build.read_spec(ROOT)
    cached=ROOT/'experiments/v5_atlas_calibrated/analysis/h1_contrasts.npz'
    cachemeta=json.loads(cached.with_suffix('.json').read_text())
    assert build.sha256(cached)==cachemeta['arrays_sha256']
    assert build.sha256(DATA/'extraction.json')==cachemeta['source_extraction_sha256']
    with np.load(cached,allow_pickle=False) as z:
        targets=z['targets'].copy(); np.testing.assert_array_equal(genes,z['genes'])
        parts=[{k:z[f'{part}_{k}'].copy() for k in ['values','errors','baselines','bulk_values','bulk_errors']} for part in [0,1]]
    batches=sorted(meta.batch.unique()); partition=meta.batch.map({b:i%2 for i,b in enumerate(batches)}).to_numpy()
    control=meta.target_gene.to_numpy()=='non-targeting'; template_rows={};reference_rows={}
    for part in [0,1]:
        templates=[];references=[]
        for batch in batches:
            rows=np.flatnonzero(control&(meta.batch.to_numpy()==batch)&(partition==part))
            if not len(rows):continue
            rng=np.random.default_rng(SEED+part+sum(batch.encode())); rows=rng.permutation(rows)
            templates.extend(rows[:len(rows)//2]);references.extend(rows[len(rows)//2:])
        template_rows[part]=np.array(templates);reference_rows[part]=np.array(references)
        assert not set(templates)&set(references)
    atlas,atlas_sha=load_atlas_model(); records=[]; profiles={}; identity=[]
    for direction in [0,1]:
        train_part=direction; test_part=1-direction
        template=x[template_rows[test_part]]; ref=x[reference_rows[test_part]]
        _,p,b,var=moments(template); ref_profile,ref_jack=profile_jackknife(ref)
        def compute(target_index, seed):
            target=str(targets[target_index]); own=int(np.flatnonzero(genes==target)[0])
            actual_rows=np.flatnonzero((meta.target_gene.to_numpy()==target)&(partition==test_part))
            n=min(400,len(actual_rows),template.shape[0])
            rng=np.random.default_rng(SEED+direction*10000+target_index*100+seed)
            tr=rng.choice(actual_rows,n,replace=False); rows=rng.choice(template.shape[0],n,replace=False)
            actual=x[tr]; truth=de(actual,ref,genes); tp,tj=profile_jackknife(actual)
            values,errors=parts[train_part]['values'][target_index],parts[train_part]['errors'][target_index]
            supported=(parts[train_part]['baselines'][target_index]>5)&(errors>0)&(genes!=target)
            specific=np.sign(values)*np.maximum(np.abs(values)-errors,0); specific[~supported]=0
            entries=[]; arrays={}
            for name,setting in CANDIDATES.items():
                effect=(.7+.8*setting['blend'])*specific
                requested=None
                if name=='v5_dual':
                    bv=parts[train_part]['bulk_values'][target_index]; be=parts[train_part]['bulk_errors'][target_index]
                    delta=np.sign(bv)*np.maximum(np.abs(bv)-be,0);delta[~supported]=0
                    requested=np.expm1(np.maximum(np.log1p(50000*b)+delta,0));requested[own]*=.1;requested/=requested.sum()
                generated,diagnostic=emit(template,rows,effect,own,SEED+direction*10000+target_index*100+seed,setting,atlas,requested)
                prediction=de(generated,ref,genes)
                metrics=direction_metrics(prediction,truth,genes,own)
                gp,gj=profile_jackknife(generated)
                _,_,_,gv=moments(generated)
                entries.append({'target_gene':target,'direction':direction,'seed':seed,'candidate':name,'cells':n,
                    'source_cells_available':len(actual_rows),**metrics,'variance_ratio_to_template':float(gv.sum()/var.sum()),
                    'mean_detected_genes':float(generated.getnnz(axis=1).mean()),**diagnostic})
                arrays[(target,direction,seed,name)]=(gp,tp,ref_profile,gj,tj,ref_jack)
            ids={'target_gene':target,'direction':direction,'seed':seed,'cells':n,
                'actual_source_ids':meta.iloc[tr].source_cell_id.tolist(),
                'template_source_ids':meta.iloc[template_rows[test_part][rows]].source_cell_id.tolist(),
                'reference_control_ids':meta.iloc[reference_rows[test_part]].source_cell_id.tolist(),
                'effect_fit_partition':train_part,'evaluation_partition':test_part}
            return entries,arrays,ids
        with ThreadPoolExecutor(max_workers=2) as pool:
            for ti in range(len(targets)):
                jobs=[pool.submit(compute,ti,seed) for seed in [0,1]]
                for job in jobs:
                    entries,arrays,ids=job.result();records.extend(entries);profiles.update(arrays);identity.append(ids)
                pd.DataFrame(records).to_csv(A/'source_generator_scores.csv',index=False)
                print(f'H1 direction {direction}: {ti+1}/{len(targets)} targets, actual rank tests complete',flush=True)
    summary=candidates_summary(records,profiles,genes,targets,official_targets)
    summary.to_csv(A/'source_candidate_summary.csv',index=False)
    # Nested outer target groups choose settings using other targets only.
    folds=target_folds(targets); outer=[]
    for fold in range(5):
        train_targets=targets[folds!=fold];test_targets=targets[folds==fold]
        train_records=[r for r in records if r['target_gene'] in set(train_targets)]
        test_records=[r for r in records if r['target_gene'] in set(test_targets)]
        selected=choose(candidates_summary(train_records,profiles,genes,train_targets,official_targets))
        test_summary=candidates_summary(test_records,profiles,genes,test_targets,official_targets)
        for row in test_summary.to_dict('records'):
            outer.append({'fold':fold,'selected_on_training_targets':row['candidate']==selected,**row})
    outer=pd.DataFrame(outer);outer.to_csv(A/'source_outer_fold_scores.csv',index=False)
    nested=outer[outer.selected_on_training_targets].mean(numeric_only=True).to_dict()
    baseline=outer[outer.candidate=='v4'].mean(numeric_only=True).to_dict()
    accepted=nested['fidelity']>=baseline['fidelity']-.005 and nested['jaccard']>=baseline['jaccard']-.005 and nested['expression_ratio']<=baseline['expression_ratio']
    selected=choose(summary) if accepted else 'v4'
    # A source-calibrated setting need not transfer; preserve V4 directions and context policy.
    setting=CANDIDATES[selected].copy()
    if setting['blend']==1.:
        # Full replacement was only a diagnostic; final response update is capped at 25%.
        setting['blend']=.25
    selection={'candidate':selected,'production_setting':setting,'nested_adoption_supported':bool(accepted),
        'nested_procedure':nested,'nested_v4_comparator':baseline,'source_targets':25,'seeds':2,'directions':2,
        'script_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/'v6_methods.py'),
        'source_sha256':source['counts_sha256'],'contrast_sha256':cachemeta['arrays_sha256'],'upstream_sha256':atlas_sha,
        'scope':'Unanchored 25-target H1 source panel, distinct controls and opposite batch-partition effects; target-held-out hyperparameters. Sizes are reduced when fewer than 400 distinct cells exist. No hidden A/B/C accuracy estimate.',
        'expression_metric_scope':'Source-panel capped-jackknife ratio; no competition anchors. Source control/template partitions differ by cell ID.',
        'response_scope':'Source gain comparison starts from expanded H1 effects in both estimators; it does not reproduce every V4/V5 context-policy/coverage difference.',
        'production_response_policy':'V4 A; B/C same-sign-only capped effect update, zero V4 coordinates retained. No neural causal counterfactual.'}
    build.write_json(A/'source_selection.json',selection);build.write_json(A/'source_cell_identity.json',identity)
    arrays={}
    for candidate in CANDIDATES:
        for d in [0,1]:
            for s in [0,1]:
                for k,key in enumerate(['pred','truth','control','pred_jack','truth_jack','control_jack']):
                    arrays[f'{candidate}_{d}_{s}_{key}']=np.stack([profiles[(str(t),d,s,candidate)][k] for t in targets])
    np.savez_compressed(A/'source_profile_cache.npz',genes=genes,targets=targets,**arrays)
    print(summary.to_string(index=False),flush=True); print('Selected source setting',selected,setting,'nested supported',accepted,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
