"""Reproduce a descriptive leaderboard gap audit; no model fitting or API calls."""
from pathlib import Path
import argparse
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/post_v29_research'
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
LABELS = ('Perturbation discrimination', 'Expression accuracy', 'DE fold-change accuracy',
          'Direction fidelity', 'Direction reach', 'Significance overlap')
EVIDENCE = ('train_v7_transition.py', 'train_v26_decoder.py', 'v20_methods.py',
            'V14_RESULTS.md', 'V18_RESULTS.md', 'V27_RESULTS.md', 'V29_RESULTS.md',
            'POST_V20_RETROSPECTIVE.md', 'PROJECT_RETROSPECTIVE.md')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build():
    snapshot_path = OUT / 'leaderboard_snapshot.json'
    ledger_path = ROOT / 'experiments/v24_cd4_response/analysis/leaderboard_result.json'
    snapshot = json.loads(snapshot_path.read_text(encoding='utf-8-sig'))
    own = json.loads(ledger_path.read_text(encoding='utf-8-sig'))
    assert tuple(own['metrics'].keys()) == METRICS
    assert abs(sum(own['metrics'].values()) / 6 - own['overall']) < 1e-12
    comparisons = []
    for team in snapshot['teams']:
        scores = team['scaled_metrics']
        assert set(scores) == set(METRICS)
        mean = sum(scores.values()) / 6
        # Screenshot components are rounded to three decimals; overall to four.
        assert abs(mean - team['overall']) <= 0.00055
        gap = {key: scores[key] - own['metrics'][key] for key in METRICS}
        gap_sum = sum(gap.values())
        assert gap_sum > 0
        comparisons.append({
            'team': team['team'], 'model': team['model'],
            'reported_overall': team['overall'], 'rounded_component_mean': mean,
            'reported_overall_gap': team['overall'] - own['overall'],
            'rounded_component_overall_gap': gap_sum / 6,
            'scaled_metric_gaps': gap,
            'arithmetic_share_of_rounded_gap': {k: gap[k] / gap_sum for k in METRICS},
            'expression_plus_lfc_share': (gap['mse'] + gap['nmae']) / gap_sum,
        })
    result = {
        'review_date': snapshot['review_date'],
        'own_entry': own['entry_id'], 'own_overall': own['overall'],
        'own_scaled_metrics': own['metrics'], 'own_raw_metrics': own['raw_metrics'],
        'own_panel': own['panel'], 'own_anchor_set': own['anchor_set'],
        'comparison_scope': 'User screenshot; competitor anchors not visible, not API verified. '
                            'Arithmetic is descriptive, not causal attribution or uncertainty.',
        'input_sha256': {'snapshot': digest(snapshot_path), 'v24_receipt': digest(ledger_path),
                         'sources_manifest': digest(OUT / 'sources.json')},
        'local_evidence_sha256': {name: digest(ROOT / name) for name in EVIDENCE},
        'comparisons': comparisons,
    }
    assert all(math.isfinite(v) for c in comparisons for v in c['scaled_metric_gaps'].values())
    lines = [
        '# Leaderboard metric gap from the supplied screenshot', '',
        'Reviewed 2026-10-15. Scaled components are higher-is-better. '
        'Competitor values are manually transcribed rounded screenshot values; '
        'their anchor IDs are not visible. These are descriptive comparisons.', '',
        '| Component | V24 | iCell | X-89 | iCell minus V24 |',
        '|---|---:|---:|---:|---:|',
    ]
    first, second = snapshot['teams'][:2]
    for key, label in zip(METRICS, LABELS):
        lines.append(f"| {label} | {own['metrics'][key]:.4f} | "
                     f"{first['scaled_metrics'][key]:.3f} | "
                     f"{second['scaled_metrics'][key]:.3f} | "
                     f"{first['scaled_metrics'][key] - own['metrics'][key]:+.4f} |")
    lines.extend(['', f"V24 overall: {own['overall']:.9f}. Screenshot iCell: "
                  f"{first['overall']:.4f}; X-89: {second['overall']:.4f}.", '',
                  f"Expression and fold-change accuracy account arithmetically for "
                  f"{comparisons[0]['expression_plus_lfc_share']:.1%} of the rounded "
                  'six-component gap to iCell. This does not identify its method or '
                  'predict the gain from changing one part of our pipeline.', '',
                  'The exact V24 raw expression error is '
                  f"{own['raw_metrics']['expr_mse_unbiased_capped_norm']:.6f}; "
                  'its raw DE fold-change error is '
                  f"{own['raw_metrics']['de_wilcoxon_lfc_nmae']:.6f}. "
                  'These raw errors are lower-is-better.', ''])
    return {'gap_summary.json': (json.dumps(result, indent=2) + '\n').encode(),
            'metric_gaps.md': '\n'.join(lines).encode()}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify deterministic archived outputs.')
    args = parser.parse_args()
    outputs = build()
    for name, body in outputs.items():
        path = OUT / name
        if args.check:
            assert path.read_bytes() == body, f'Outdated research artifact: {path}'
        else:
            path.write_bytes(body)
    print('Leaderboard arithmetic and evidence checks passed.' if args.check
          else 'Wrote descriptive leaderboard comparison; no predictions changed.')

if __name__ == '__main__':
    main()
