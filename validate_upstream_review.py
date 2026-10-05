"""Validate upstream-review provenance, folds, coverage tables and unchanged V4 model."""
from pathlib import Path
import gzip
import importlib.util
import json
import re
import sys
import numpy as np
import pandas as pd
import build_submission as build
from review_upstream_repositories import git
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'analysis/upstream_review'


def main():
    review=json.loads((OUT/'repository_review.json').read_text())
    assert review['script_sha256']==build.sha256(ROOT/'review_upstream_repositories.py')
    for repo in review['repositories']:
        assert git(repo['repository'],'rev-parse','HEAD').strip()==repo['head']
        assert not git(repo['repository'],'status','--porcelain').strip()
        for record in repo['tracked_files']:assert build.sha256(ROOT/repo['repository']/record['path'])==record['sha256']
    centering=json.loads((OUT/'centering_validation.json').read_text());seen=[];records=centering['nested_folds']
    for record in records:
        assert not set(record['train_targets'])&set(record['test_targets']);seen.extend(record['test_targets'])
    assert len(seen)==len(set(seen))==centering['targets']
    all_targets=set(seen)
    for record in records:
        assert set(record['train_targets'])|set(record['test_targets'])==all_targets
    n=sum(r['contrasts'] for r in records);assert n==centering['contrasts']
    outer=pd.read_csv(OUT/'centering_outer_scores.csv');candidates=pd.read_csv(OUT/'centering_candidates.csv')
    assert len(outer)==5*len(centering['grid']) and not outer.duplicated(['fold','candidate']).any()
    assert len(candidates)==len(centering['grid'])
    for candidate in candidates.to_dict('records'):
        group=outer[outer.candidate==candidate['candidate']]
        assert set(group.fold)==set(range(5)) and int(group.contrasts.sum())==n
        assert json.loads(candidate['setting'])==centering['grid'][candidate['candidate']]
        for field in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy']:
            np.testing.assert_allclose(candidate[field],np.average(group[field],weights=group.contrasts),atol=1e-12)
    baseline=candidates[candidates.setting.map(lambda value:json.loads(value)=={'gain':1.,'centering':0.})].iloc[0]
    for field in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy']:
        np.testing.assert_allclose(baseline[field],centering['v4_fixed_baseline'][field],atol=1e-12)
    for field in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy']:
        np.testing.assert_allclose(centering['nested_procedure'][field],sum(r[field]*r['contrasts'] for r in records)/n,atol=1e-12)
    np.testing.assert_allclose(centering['nested_procedure']['relative_zero_loss'],centering['nested_procedure']['loss']/centering['nested_procedure']['zero_loss'],atol=1e-12)
    assert centering['script_sha256']==build.sha256(ROOT/'evaluate_atlas_transfer.py')
    probe=json.loads((OUT/'h1_coverage_probe.json').read_text());coverage=pd.read_csv(OUT/'h1_public_split_coverage.csv')
    assert len(coverage)==900 and not coverage.duplicated(['split','target_gene']).any()
    combined=coverage.pivot(index='target_gene',columns='split',values='available_cells').fillna(0)
    assert int((combined.sum(axis=1)>0).sum())==probe['all_splits_overlap']
    assert int((combined.train>0).sum())==probe['training_overlap']
    assert set(combined.index[(combined.train==0)&(combined.sum(axis=1)>0)])==set(probe['additional_targets'])
    expected=combined.astype(int).copy();expected['all_splits_cells']=expected.sum(axis=1)
    expected['new_beyond_training']=(expected.train==0)&(expected.all_splits_cells>0)
    actual=pd.read_csv(OUT/'h1_public_target_coverage.csv',index_col='target_gene');actual.columns.name='split'
    pd.testing.assert_frame_equal(actual,expected,check_dtype=False)
    assert probe['script_sha256']==build.sha256(ROOT/'probe_atlas_sources.py')
    for source in probe['sources']:
        cache=ROOT/'reference/control_profiles/ranges_H1' if source['split']=='train' else ROOT/'reference/upstream_review'/('ranges_H1_'+source['split'])
        assert json.loads((cache/'identity.json').read_text())==source['source']
    model=json.loads((ROOT/'experiments/v4_lineage_response/analysis/reference_model.json').read_text())
    assert build.sha256(ROOT/'experiments/v4_lineage_response/analysis/reference_model_arrays.npz')==model['model_arrays_sha256']
    assert build.sha256(ROOT/'fit_lineage_submission.py')==model['script_sha256']
    history=json.loads(git('Virtual-Cell-Challenge-2026','show','8257ada:docs/v38-inputs.sha256.json'))['files']
    original={r['file']:r['sha256'] for r in model['context_original_inputs']}
    original.update({name:build.sha256(ROOT/name) for name in ['gene_names.csv','pert_counts.csv']})
    matches={name:{'our_sha256':value,'historical_sha256':history[name],'match':value==history[name]} for name,value in original.items()}
    assert all(r['match'] for r in matches.values())
    promoter=json.loads((OUT/'promoter_probe.json').read_text());assert promoter['script_sha256']==build.sha256(ROOT/'probe_atlas_sources.py')
    annotation=ROOT/'reference/upstream_review'/promoter['annotation']['filename'];assert build.sha256(annotation)==promoter['annotation']['sha256']
    table=[]
    with gzip.open(annotation,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if fields[2]!='gene':continue
            attributes=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            if 'gene_name' in attributes:table.append({'gene':attributes['gene_name'],'chromosome':fields[0],
                'strand':fields[6],'tss':int(fields[3] if fields[6]=='+' else fields[4])})
    spec=importlib.util.spec_from_file_location('reviewed_atlas_prepare',ROOT/'Virtual-Cell-Challenge-2026/prepare.py')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    pairs=module.pairs(pd.DataFrame(table),pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy(),pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy())
    ours=pd.read_csv(OUT/'promoter_neighbor_pairs.csv')
    pd.testing.assert_frame_equal(ours,pairs[ours.columns],check_dtype=False)
    import hashlib
    reproduced_sha=hashlib.sha256(pairs.to_csv(index=False,lineterminator='\n').encode()).hexdigest()
    promoter_identity={'our_logical_pairs_match_reviewed_upstream':True,'reproduced_csv_sha256':reproduced_sha,
        'historical_csv_sha256':history['official_pairs.csv'],'historical_csv_byte_match':reproduced_sha==history['official_pairs.csv']}
    observations=pd.read_csv(OUT/'cd4_promoter_neighbor_observations.csv');significant=observations.q_value<.05
    assert int((significant&(observations.log2fc<0)).sum())==promoter['significant_negative']
    assert int((significant&(observations.log2fc>0)).sum())==promoter['significant_positive']
    moments=pd.read_csv(OUT/'control_moment_diagnostics.csv');assert len(moments)==9 and (moments.cells==400).all() and moments.positive_integer_counts.all()
    for field in moments.select_dtypes(include='number').columns:
        assert np.isfinite(moments[field].dropna()).all()
    assert json.loads((OUT/'control_moment_scope.json').read_text())['script_sha256']==build.sha256(ROOT/'evaluate_atlas_transfer.py')
    assert 'Ran 4 tests' in (OUT/'test_results.txt').read_text(encoding='utf-8-sig') and '\nOK' in (OUT/'test_results.txt').read_text(encoding='utf-8-sig')
    build.write_json(OUT/'historical_input_identity.json',{'competition_inputs':matches,'promoter_annotation':promoter_identity,
        'scope':'Control hashes compared with previously verified V4 input hashes; gene/target files freshly hashed. Does not establish equality of hidden references or leaderboard anchor versions.'})
    package=ROOT/'experiments/v4_lineage_response/prediction.vcc';ready=json.loads((ROOT/'experiments/v4_lineage_response/analysis/submission_ready.json').read_text())
    assert package.stat().st_size==ready['bytes']
    artifact_names=['repository_review.json','centering_validation.json','control_moment_diagnostics.csv','control_moment_scope.json',
        'h1_coverage_probe.json','h1_public_split_coverage.csv','h1_public_target_coverage.csv','promoter_probe.json',
        'promoter_neighbor_pairs.csv','cd4_promoter_neighbor_observations.csv','historical_input_identity.json','test_results.txt','centering_candidates.csv','centering_outer_scores.csv',
        'Virtual-Cell-Challenge-2026_commits.csv','arc-virtual-cell-atlas_commits.csv']
    result={'passed':True,'nested_targets_checked':len(seen),'coverage_rows_checked':len(coverage),'moment_groups_checked':len(moments),
        'upstream_promoter_pairs_independently_matched':len(pairs),'v4_effect_cache_sha256_unchanged':model['model_arrays_sha256'],
        'v4_package_size_unchanged':True,'artifact_sha256':{name:build.sha256(OUT/name) for name in artifact_names},
        'root_code_sha256':{name:build.sha256(ROOT/name) for name in ['review_upstream_repositories.py','probe_atlas_sources.py','evaluate_atlas_transfer.py','test_upstream_review.py','validate_upstream_review.py']},
        'scope':'Provenance, stored numerical summaries, disjoint folds, metadata coverage, full original GTF pair reconstruction and focused tests. Does not independently recompute hidden scores or reread all source expression counts.'}
    build.write_json(OUT/'validation.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
