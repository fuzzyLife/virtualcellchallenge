"""Record the user-reported V4 publication, verify it via the read-only CLI, and compare scores."""
from datetime import datetime, timezone
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / 'experiments/v4_lineage_response/analysis'
ENTRY = 'w5criCkgkCWj1QOgP08E'
REPORTED = {
    'overall': -0.05067391016092112,
    'pds': 0.25530201756101617, 'mse': 0.0,
    'nmae': 0.06753871251359793, 'fid': -0.6838546480756342,
    'reach': 0.08684858992027715, 'jac': -0.02987813288478362,
}
RAW = {
    'pds_cosine': True, 'expr_mse_unbiased_capped_norm': False,
    'de_wilcoxon_lfc_nmae': False, 'de_wilcoxon_direction_fidelity_yield_raw': True,
    'de_wilcoxon_direction_reach_raw': True, 'de_wilcoxon_sig_jaccard': True,
}


def main():
    process = subprocess.run([sys.executable, '-m', 'vcc', '--json', 'status', ENTRY],
                             capture_output=True, text=True, encoding='utf-8', check=True, timeout=60)
    api = json.loads(process.stdout)
    assert api['entry_id'] == ENTRY and api['status'] == 'published' and api['is_terminal']
    for name, reported in REPORTED.items():
        np.testing.assert_allclose(api['score_avg' if name == 'overall' else 'score_' + name],
                                   reported, rtol=0, atol=1e-12)
    first = json.loads((ROOT / 'analysis/leaderboard_result.json').read_text())
    previous = json.loads((ROOT / 'experiments/v2_cd4_hybrid/analysis/leaderboard_result.json').read_text())
    previous_api = json.loads((ROOT / 'experiments/v2_cd4_hybrid/analysis/published_api_status.json').read_text())
    first_api = json.loads((ROOT / 'analysis/model_audit/published_api_status.json').read_text())
    ready = json.loads((ANALYSIS / 'submission_ready.json').read_text())
    assert api['partition'] == previous['partition'] == first['partition']
    assert api['panel_id'] == previous['panel'] == first['panel']
    assert api['anchor_version'] == previous['anchor_set'] == first['anchor_set']
    comparison = []
    for metric, value in REPORTED.items():
        old = previous['overall'] if metric == 'overall' else previous['metrics'][metric]
        original = first['overall'] if metric == 'overall' else first['metrics'][metric]
        comparison.append({'metric': metric, 'first': original, 'v2': old, 'v4': value,
                           'change_from_v2': value - old, 'change_from_first': value - original,
                           'overall_contribution_change_from_v2': None if metric == 'overall' else (value - old) / 6})
    comparison = pd.DataFrame(comparison)
    np.testing.assert_allclose(comparison.iloc[1:].overall_contribution_change_from_v2.sum(),
                               REPORTED['overall'] - previous['overall'], rtol=0, atol=1e-12)
    raw = pd.DataFrame([{'metric': name, 'first': first_api[name], 'v2': previous_api[name],
                        'v4': api[name], 'change_from_v2': api[name] - previous_api[name],
                        'higher_is_better': higher,
                        'improved_from_v2': api[name] > previous_api[name] if higher else api[name] < previous_api[name]}
                       for name, higher in RAW.items()])
    expression_change = api['expr_mse_unbiased_capped_norm'] / previous_api['expr_mse_unbiased_capped_norm'] - 1
    gain = REPORTED['overall'] - previous['overall']
    fidelity_share = (REPORTED['fid'] - previous['metrics']['fid']) / 6 / gain
    result = {
        'source': 'User-provided CLI publication output, independently confirmed by read-only vcc --json status.',
        'reported_date': '2026-10-04', 'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'entry_id': ENTRY, 'model_name': api['model_name'], 'final_status': 'published',
        'rank_at_publication': 977, 'rank_at_api_check': api['rank'],
        'overall': REPORTED['overall'], 'metrics': {k: v for k, v in REPORTED.items() if k != 'overall'},
        'raw_metrics': {k: api[k] for k in RAW}, 'partition': api['partition'], 'panel': api['panel_id'],
        'anchor_set': api['anchor_version'], 'package_sha256': ready['sha256'],
        'upload_bytes': ready['bytes'], 'upload_duration_seconds': 800,
        'overall_change_from_v2': gain, 'overall_change_from_first': REPORTED['overall'] - first['overall'],
        'raw_expression_error_relative_change_from_v2': expression_change,
        'fidelity_fraction_of_overall_gain_from_v2': fidelity_share,
        'interpretation': 'All five non-expression scaled metrics improve; raw expression error worsens while its scaled score stays at zero. Overall remains negative; fidelity is the largest negative component. Multiple model changes prevent attributing the gain to one change. V3 has no recorded leaderboard score.',
        'readiness_record_scope': 'submission_ready.json is the immutable generation-time snapshot; this separate publication record is the current upload/scoring evidence.',
    }
    (ANALYSIS / 'published_api_status.json').write_text(json.dumps(api, indent=2) + '\n', encoding='utf-8')
    (ANALYSIS / 'leaderboard_result.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    comparison.to_csv(ANALYSIS / 'leaderboard_comparison.csv', index=False)
    raw.to_csv(ANALYSIS / 'raw_metric_comparison.csv', index=False)
    from audit_v4_strategy import document_results
    document_results()
    print(comparison.to_string(index=False))
    print(f'Raw expression error increased {expression_change:.2%}; fidelity accounts for {fidelity_share:.2%} of the aggregate gain.')
    print('Recorded publication rank 977; current API rank', api['rank'])


if __name__ == '__main__':
    main()
