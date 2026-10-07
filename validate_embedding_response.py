"""Independently check fold boundaries, training masks and cell evaluation provenance."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import build_submission as build
from train_embedding_response import OUT, ROOT
from embedding_response_methods import VARIANTS


def main():
    source=json.loads((OUT/'source_provenance.json').read_text())
    fit=json.loads((OUT/'fit_manifest.json').read_text())
    for name,digest in source['sources'].items():assert build.sha256(ROOT/name)==digest,name
    assert build.sha256(OUT/'source_cache.npz')==source['cache_sha256']
    assert build.sha256(ROOT/'train_embedding_response.py')==fit['script_sha256']==source['script_sha256']
    assert build.sha256(ROOT/'embedding_response_methods.py')==fit['methods_sha256']
    assert build.sha256(OUT/'h1_lineage_predictions.npz')==fit['prediction_cache_sha256']
    with np.load(OUT/'source_cache.npz',allow_pickle=False) as z:
        genes,targets=z['genes'],z['targets'];v=z['valid'];families=z['row_families'];row_targets=z['row_targets'];n=len(row_targets)
        assert genes.shape==(18533,) and targets.shape==(300,)
        assert not v[:,np.isin(genes,targets)].any()
        assert np.isfinite(z['values']).all() and (z['values'][~v]==0).all()
        assert np.isfinite(z['errors']).all() and (z['errors'][v]>0).all()
        assert not v[z['row_contexts']==4][:, z['profiles'][4]<=5].any()
        assert len(set(genes))==len(genes)
    splits=json.loads((OUT/'split_manifest.json').read_text());assert len(splits)==6
    for s in splits:
        a,b=set(s['train_rows']),set(s['test_rows']);assert not a&b and a|b==set(range(n))
        if s['regime']=='lineage_holdout':assert not set(s['train_families'])&set(s['test_families'])
        else:assert not set(s['train_targets'])&set(s['test_targets'])
    score=pd.read_csv(OUT/'contrast_scores.csv')
    assert not score.duplicated(['regime','context','target_gene','variant']).any()
    assert set(score.variant)==set(VARIANTS) and len(score)==n*len(VARIANTS)*2
    assert (score.genes_evaluated>0).all()
    decision_path=OUT/'adoption_decision.json'
    if decision_path.exists():
        decision=json.loads(decision_path.read_text())
        assert build.sha256(ROOT/'evaluate_embedding_response.py')==decision['count_script_sha256']
        assert build.sha256(ROOT/'v6_methods.py')==decision['vcc_methods_sha256']
        assert build.sha256(ROOT/'embedding_response_methods.py')==decision['head_methods_sha256']
        count=pd.read_csv(OUT/'generated_cell_scores.csv')
        expected=len(decision['primary_targets'])*decision['seeds']*decision['partitions']*5
        assert len(count)==expected
        assert not count.duplicated(['target_gene','partition','seed','candidate']).any()
        assert count.groupby(['target_gene','partition','seed']).cells.nunique().max()==1
        assert (count.cells>=200).all() and np.isfinite(count.fidelity).all()
        ids=json.loads((OUT/'generated_cell_identity.json').read_text())
        assert len(ids)==expected//5
        for item in ids:
            a,b,c=map(set,[item['actual_source_ids'],item['template_source_ids'],item['reference_control_ids']])
            assert len(a)==len(b)==item['n'] and not a&b and not a&c and not b&c
        screen=pd.read_csv(OUT/'contrast_summary.csv')
        for candidate, guards in decision['decisions'].items():
            for family, expected in guards['summary_context_guards'].items():
                subset=screen[(screen.regime=='lineage_holdout')&(screen.family==family)]
                a=subset[subset.variant==candidate].iloc[0];b=subset[subset.variant=='identity_rank'].iloc[0]
                actual=bool(a.mse<=b.mse+1e-12 and a.fidelity_proxy>=b.fidelity_proxy-1e-12 and a.jaccard_proxy>=b.jaccard_proxy-1e-12)
                assert actual==expected, (candidate,family)
        proof_path=OUT/'mask_boundary_verification.json'
        if proof_path.exists():
            proof=json.loads(proof_path.read_text())
            assert proof['all_arrays_exactly_equal']
            assert proof['old_cache_sha256']==proof['current_cache_sha256']==fit['prediction_cache_sha256']
        if decision['selected']=='v4':assert not decision['package_generated']
    report={'passed':True,'source_rows':n,'contrast_records':len(score),'holdout_splits':len(splits),
        'count_evaluation_complete':decision_path.exists(),'fit_manifest_sha256':build.sha256(OUT/'fit_manifest.json'),
        'validator_sha256':build.sha256(Path(__file__))}
    build.write_json(OUT/'validation.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
