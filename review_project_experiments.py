"""Audit historical evidence without fitting models or changing submissions.

Writes only small tables under analysis/project_retrospective/.
PROJECT_RETROSPECTIVE.md contains the human-reviewed interpretation.
"""
from __future__ import annotations
import hashlib
import json
import platform
import re
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/project_retrospective'
INPUTS = {}
METRICS = ['pds', 'mse', 'nmae', 'fid', 'reach', 'jac']
RAW = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae',
       'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw',
       'de_wilcoxon_sig_jaccard']
# Curated scope, not an inference that every plan in a document was executed.
CATALOG = [
 ('Controls and rank analyses', 'Full/300-gene clustering, external mapping and marker audit', 'Baseline separation; no causal response labels', 'diagnostic', 'TARGET_LINEAGE_INTERPRETATION.md;EXTERNAL_RANK_MAPPING_RESULTS.md;TARGET_CELL_DISPERSION_RESULTS.md'),
 ('V1/V2', 'K562 transfer, then donor-calibrated CD4 and zero-aware count increases', 'Published improvement, still one shared response', 'published', 'COMPETITION_ANALYSIS.md;V2_SUBMISSION.md;build_submission.py'),
 ('V3', 'Basal NNLS ensembles, CD4/H1 fitting and common templates', 'Generated and audited; no recorded score', 'unsubmitted', 'V3_SUBMISSION.md'),
 ('V4', 'Lineage-guided source policies and downstream response coverage', 'Large published gain; raw expression worsened; bundled changes', 'published', 'V4_MODEL_STRATEGY.md;V4_SUBMISSION.md'),
 ('V5/V9', 'Mean/bulk moment calibration; activation-preserving partial calibration', 'Five scaled gains versus V4 outweighed by fidelity loss', 'published', 'V5_MODEL_STRATEGY.md;V5_LEADERBOARD_RESULTS.md;V9_LEADERBOARD_RESULTS.md'),
 ('Control reconstruction', '36 ridge candidates; 24 shallow/deep/bottleneck/residual neural networks', 'Neural raw MSE improves; causal intervention accuracy untested', 'observational only', 'GENE300_MODEL_RESULTS.md;GENE300_NEURAL_RESULTS.md;GENE300_NEURAL_MODELS.md'),
 ('V6', 'Response calibration, fidelity guards and neural-template gate', 'Published below V4; local guards did not protect the leaderboard', 'published', 'V6_MODEL_STRATEGY.md;V6_SUBMISSION.md;V6_LEADERBOARD_RESULTS.md'),
 ('Embedding pilot', 'Four Tahoe/PCA response heads and H1 generated-cell blends', 'Tested blends failed; proteomics checkpoint absent', 'rejected limited heads', 'CROSS_MODAL_MODEL_REVIEW.md;EMBEDDING_RESPONSE_PILOT.md'),
 ('V7', 'Expanded HCT/HEK, nested pooling, supervised neural/ridge alternatives', 'Alternatives rejected; upload equals V4 effects with new sampling', 'published fallback', 'V7_MODEL_STRATEGY.md;V7_LEADERBOARD_RESULTS.md'),
 ('V8', 'Common-origin PDS repair, response and count calibration', 'Diagnostic package; no recorded score', 'unsubmitted', 'TRANSFER_VALIDATION_AUDIT.md;V8_SUBMISSION.md'),
 ('V10', 'Measured family pooling, dose correction, EB shrinkage and compact emission', 'Major published gain; response and emission changed together', 'published', 'V10_STRATEGY.md;V10_LEADERBOARD_RESULTS.md;v10_methods.py'),
 ('V11', 'Consensus rescue and noise 0.65, isolated/combined local ablations', 'Published gain; expression worsened; population bias remains', 'published', 'V11_STRATEGY.md;V11_LEADERBOARD_RESULTS.md;v11_methods.py'),
 ('V12 distributions', 'Empirical/gamma mixtures, response-linked variation, thinning and higher Poisson noise', 'Better null realism often loses intervention fidelity/error guards', 'rejected', 'V12_RESEARCH.md;v12_methods.py'),
 ('V12 sources and priors', 'K562 Flex, six platform variants, heterogeneity/amplitude and Reactome priors', 'No accepted candidate; pathway screen did not beat shuffled annotations', 'rejected limited priors', 'V12_RESEARCH.md;v12_harmonized_methods.py;screen_v12_responses.py'),
 ('V13 programs', 'Low-rank response bases: ranks 8/32/64 and blends', 'Count finalists failed; not the production response', 'rejected', 'V13_STRATEGY.md;v13_methods.py'),
 ('V13 robust', 'Median-radius clipping and Student-t source reliability', 'Student3 improves five primary development metrics; expression guard fails', 'unconfirmed lead', 'V13_STRATEGY.md;v13_robust_methods.py;assess_v13_submission.py'),
 ('V13 quality', 'Unknown-KD correction, soft weights and hard KD gate', 'Unknown-KD released; incomplete confirmation; tiny published gain', 'published exploratory', 'V13_STRATEGY.md;V13_LEADERBOARD_RESULTS.md;v13_quality_methods.py'),
 ('Data feasibility', 'Reactome/fingerprinting code, RPE1 overlap, recovered construct KD', 'Fingerprint fitting not run; RPE1 zero target overlap; actual KD unused after freeze', 'mixed / partly untested', 'POST_V11_RESEARCH.md;audit_v13_rpe1.py;audit_v13_knockdown.py'),
]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def track(relative):
    p = ROOT / relative
    key = p.relative_to(ROOT).as_posix()
    digest = sha(p)
    if key in INPUTS and INPUTS[key] != digest:
        raise ValueError('Input changed while reading: ' + key)
    INPUTS[key] = digest
    return p


def read_json(relative):
    return json.loads(track(relative).read_text(encoding='utf-8-sig'))


def read_csv(relative):
    return pd.read_csv(track(relative))


def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def export(name, rows):
    pd.DataFrame(rows).to_csv(OUT / name, index=False, lineterminator='\n')


def published():
    files = [ROOT / 'analysis/leaderboard_result.json', *ROOT.glob('experiments/*/analysis/leaderboard_result.json')]
    rows = []
    for p in files:
        d = read_json(p.relative_to(ROOT))
        version = 1 if p.parent == ROOT / 'analysis' else int(re.match(r'v(\d+)_', p.parent.parent.name)[1])
        assert d['final_status'] == 'published'
        assert abs(sum(d['metrics'][m] for m in METRICS) / 6 - d['overall']) < 1e-12
        rows.append({'version': version, 'label': 'First' if version == 1 else f'V{version}',
                     'entry_id': d['entry_id'], 'overall': d['overall'], **d['metrics'],
                     **{k: d.get('raw_metrics', {}).get(k) for k in RAW},
                     'panel': d['panel'], 'partition': d['partition'], 'anchor_set': d['anchor_set'],
                     'record': p.relative_to(ROOT).as_posix()})
    rows.sort(key=lambda r: r['version'])
    assert len(rows) == 10 and len({r['entry_id'] for r in rows}) == 10
    assert len({(r['panel'], r['partition'], r['anchor_set']) for r in rows}) == 1
    best = -float('inf')
    for i, r in enumerate(rows):
        r['change_from_previous_upload'] = None if i == 0 else r['overall'] - rows[i - 1]['overall']
        r['change_from_previous_best'] = None if i == 0 else r['overall'] - best
        best = max(best, r['overall'])
        r['best_so_far'] = best
    export('published_history.csv', rows)
    by = {r['version']: r for r in rows}
    changes = []
    for a, b in [(1, 2), (2, 4), (4, 5), (4, 6), (4, 7), (4, 9), (4, 10), (10, 11), (11, 13)]:
        parts = [(by[b][m] - by[a][m]) / 6 for m in METRICS]
        assert abs(sum(parts) - (by[b]['overall'] - by[a]['overall'])) < 1e-12
        for m, value in zip(METRICS, parts):
            changes.append({'from': by[a]['label'], 'to': by[b]['label'], 'component': m,
                            'scaled_change': value * 6, 'contribution_to_overall_change': value})
    export('component_attribution.csv', changes)
    return rows


def geometry():
    folder = Path('experiments/v13_response_programs/analysis')
    model = read_json(folder / 'reference_model.json')
    p = track(folder / 'reference_model_arrays.npz')
    assert sha(p) == model['model_arrays_sha256']
    for name, expected in model['dependencies'].items():
        assert sha(track(name)) == expected, 'Changed frozen dependency: ' + name
    with np.load(p, allow_pickle=False) as z:
        genes, targets, effects = z['genes'], z['targets'], z['effects']
    assert effects.shape == (3, 300, 18533) and np.isfinite(effects).all()
    profiles = []
    for context in 'ABC':
        f = read_csv(folder / f'gene_statistics_{context}.csv')
        np.testing.assert_array_equal(f.gene_name.to_numpy(), genes)
        profiles.append(f.mean_cpm.to_numpy())
    q = np.stack(profiles)
    own = np.array([np.flatnonzero(genes == t)[0] for t in targets])
    gate = np.clip(q / 5, 0, 1)
    all_gene = (q >= 5).all(axis=0)
    all_target = (q[:, own] >= 5).all(axis=0)
    pairs = []
    for a, b in combinations(range(3), 2):
        both = (effects[a] != 0) & (effects[b] != 0)
        reversals = both & (np.signbit(effects[a]) != np.signbit(effects[b]))
        diff = np.abs(effects[a][all_target][:, all_gene] - effects[b][all_target][:, all_gene])
        row = {'contexts': 'ABC'[a] + '-' + 'ABC'[b], 'joint_nonzero_effects': int(both.sum()),
               'opposite_sign_effects': int(reversals.sum()), 'exact_equal_fraction': float((effects[a] == effects[b]).mean()),
               'max_difference_fully_open_gates': float(diff.max())}
        assert row['opposite_sign_effects'] == 0 and row['max_difference_fully_open_gates'] == 0
        pairs.append(row)
    target_rows = []
    for ti, target in enumerate(targets):
        scales = gate * gate[:, own[ti]][:, None]
        active = (scales > 0).all(axis=0)
        recovered = effects[:, ti, active] / scales[:, active]
        maximum = float(np.ptp(recovered, axis=0).max()) if active.any() else 0.
        assert maximum < 2e-6
        target_rows.append({'target': target, 'target_gate_A': float(gate[0, own[ti]]),
                            'target_gate_B': float(gate[1, own[ti]]), 'target_gate_C': float(gate[2, own[ti]]),
                            'max_gate_removed_disagreement': maximum})
    export('context_effect_pairs.csv', pairs)
    export('target_context_gates.csv', target_rows)
    result = {'fully_open_output_genes': int(all_gene.sum()), 'fully_open_targets': int(all_target.sum()),
              'identical_downstream_target_gene_pairs': int(all_gene.sum() * all_target.sum() - all_target.sum()),
              'max_gate_removed_disagreement': max(r['max_gate_removed_disagreement'] for r in target_rows),
              'pairs': pairs, 'scope': 'Pre-emission log2 effects; final normalized fold changes/DE calls can differ. No hidden truth comparison.'}
    write_json('context_effect_geometry.json', result)
    return result


def local_candidates():
    rows = []
    stages = [('v11', 'experiments/v11_consensus_coverage/analysis'),
              ('v13_robust', 'experiments/v13_response_programs/analysis/robust_stage'),
              ('v13_quality', 'experiments/v13_response_programs/analysis/quality_stage')]
    for stage, folder in stages:
        frame = read_csv(folder + '/development_summary.csv')
        rows.extend({'stage': stage, **r} for r in frame.to_dict('records'))
    export('selected_local_ablation_history.csv', rows)
    robust = pd.DataFrame([r for r in rows if r['stage'] == 'v13_robust' and r['reference_mode'] == 'full_independent']).set_index('candidate')
    fields = ['pds', 'fidelity', 'reach', 'jaccard', 'nmae', 'common_expression_ratio']
    delta = robust.loc['student3', fields] - robust.loc['v11', fields]
    assert sum(delta[k] > 0 if k in fields[:4] else delta[k] < 0 for k in fields) == 5
    anchor = read_csv('experiments/v13_response_programs/analysis/robust_stage/development_source_anchored_summary.csv').set_index('candidate')
    result = {'primary_metric_changes': {k: float(delta[k]) for k in fields},
              'relative_expression_error_increase': float(delta.common_expression_ratio / robust.loc['v11', 'common_expression_ratio']),
              'source_anchored_development_overall_v11': float(anchor.loc['v11', 'overall']),
              'source_anchored_development_overall_student3': float(anchor.loc['student3', 'overall']),
              'scope': 'Reused development; first-seed public anchors with small denominators. No independent count confirmation or private leaderboard result for Student3.'}
    assessment = read_json('experiments/v13_response_programs/analysis/robust_stage/development_assessment.json')
    failures = {mode: [key for key, passed in record['checks'].items() if not passed]
                for mode, record in assessment['results']['student3']['modes'].items()}
    assert all(items == ['expression_improves'] for items in failures.values())
    result['failed_frozen_guards'] = failures
    write_json('student3_policy_tradeoff.json', result)
    a = read_json('experiments/v13_response_programs/analysis/quality_stage/validation_plan.json')
    b = read_json('experiments/v12_distribution_calibrated/analysis/harmonized_stage/validation_plan.json')
    def cases(plan, split):
        return {(c['family'], c['context'], x['target'], x['seed']) for c in plan['contexts'] for x in c['cells'] if x['split'] == split}
    ca, cb = cases(a, 'confirmation'), cases(b, 'confirmation')
    assert ca <= cb
    export('supplemental_confirmation_identities.csv', [dict(zip(['family', 'context', 'target', 'seed'], x)) for x in sorted(cb - ca)])
    power = {'v13_development_target_seed_cases': len(cases(a, 'development')),
             'v13_confirmation_target_seed_cases': len(ca), 'v13_confirmation_families': sorted({x[0] for x in ca}),
             'v12_supplemental_confirmation_cases_not_in_v13_plan': len(cb - ca),
             'scope': 'Identity comparison only, not a claim of unseen biology or adequate DE power. Original V13 confirmation has been used.'}
    write_json('validation_scope.json', power)
    read_json('experiments/v13_response_programs/analysis/confirmation_result.json')
    read_json('experiments/v13_response_programs/analysis/published_knockdown_audit.json')
    return result, power


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    track(Path(__file__).name)
    track('PROJECT_RETROSPECTIVE.md')
    inventory = []
    for area, scope, outcome, status, evidence in CATALOG:
        for name in evidence.split(';'):
            track(name)
        inventory.append({'area': area, 'tested_scope': scope, 'outcome': outcome, 'evidence_status': status, 'evidence': evidence})
    export('experiment_inventory.csv', inventory)
    for path in ['analysis/post_v11_research/transfer_family_summary.csv', 'analysis/post_v11_research/uncertainty_bins.csv',
                 'analysis/post_v11_research/neutral_summary.csv', 'EVALUATION_CRITERIA_AUDIT.md', 'FAQ_AUDIT.md',
                 'experiments/v12_distribution_calibrated/analysis/neutral/neutral_summary.csv']:
        track(path)
    rows = published()
    geo = geometry()
    student, power = local_candidates()
    by = {r['version']: r for r in rows}
    write_json('key_findings.json', {
        'best_version': 13, 'best_overall': by[13]['overall'],
        'v4_to_v10_gain': by[10]['overall'] - by[4]['overall'],
        'v10_to_v11_gain': by[11]['overall'] - by[10]['overall'],
        'v11_to_v13_gain': by[13]['overall'] - by[11]['overall'],
        'v13_pds_fraction_of_net_overall': by[13]['pds'] / 6 / by[13]['overall'],
        'v13_scaled_component_contributions': {k: by[13][k] / 6 for k in METRICS},
        'raw_metric_scope': 'Raw fidelity is budget-normalized yield, not sign precision; Jaccard is not recall.',
        'context_geometry': geo, 'student3': student, 'validation_scope': power})
    for name, expected in INPUTS.items():
        if sha(ROOT / name) != expected:
            raise ValueError('Input mutated during retrospective: ' + name)
    outputs = {p.relative_to(ROOT).as_posix(): sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'audit_manifest.json'}
    write_json('audit_manifest.json', {
        'scope': 'Retrospective review; no fitting, selection, evaluator reruns or upload.',
        'checks_passed': ['ten unique published entries with shared panel/anchors', 'overall equals six scaled component mean',
                          'metric contributions sum to observed changes', 'V13 model array and frozen dependencies match',
                          'all pairwise effect sign/gate geometry verified', 'Student3 five-metric tradeoff recomputed',
                          'validation-plan identities compared', 'all read input hashes unchanged'],
        'inputs': INPUTS, 'outputs': outputs,
        'runtime': {'python': platform.python_version(), 'executable': sys.executable, 'numpy': np.__version__,
                    'numpy_path': np.__file__, 'pandas': pd.__version__, 'pandas_path': pd.__file__}})
    print(json.dumps({'published_submissions': len(rows), 'verified_inputs': len(INPUTS),
                      'context_geometry': geo, 'student3': student, 'validation_scope': power}, indent=2))


if __name__ == '__main__':
    main()
