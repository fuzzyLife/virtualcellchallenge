"""Record V5 publication with read-only API confirmation and reproducible comparisons."""
from datetime import datetime, timezone
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / 'experiments/v5_atlas_calibrated/analysis'
ENTRY = 'BPI0n52HTsG0EYTFOf9H'
REPORTED = {'overall': -0.1136367314348, 'pds': 0.4478561566429275,
    'mse': 0.03529351682041667, 'nmae': 0.07018875085396455,
    'fid': -1.3291001637658362, 'reach': 0.12106109775869418,
    'jac': -0.027119746918966695}
RAW = {'pds_cosine': True, 'expr_mse_unbiased_capped_norm': False,
    'de_wilcoxon_lfc_nmae': False, 'de_wilcoxon_direction_fidelity_yield_raw': True,
    'de_wilcoxon_direction_reach_raw': True, 'de_wilcoxon_sig_jaccard': True}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def document_publication():
    """Update the publication layer; callable after regeneration of the V5 methods document."""
    record = json.loads((ANALYSIS / 'leaderboard_result.json').read_text(encoding='utf-8'))
    scaled = pd.read_csv(ANALYSIS / 'leaderboard_comparison.csv')
    raw = pd.read_csv(ANALYSIS / 'raw_metric_comparison.csv')
    lines = ['# V5 published leaderboard result', '',
        f"Entry `{ENTRY}` is **published**, with overall **{record['overall']:.6f}**. V4 remains the best project submission at **-0.050674**. "
        'All four published entries use the same validation panel and anchor set; V3 has no recorded score.', '',
        '| Scaled metric | First | V2 | V4 | V5 | V5 minus V4 | Contribution to overall change |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for row in scaled.itertuples():
        contribution = '' if row.metric == 'overall' else f'{row.overall_contribution_change_from_v4:+.6f}'
        lines.append(f'| {row.metric} | {row.first:.6f} | {row.v2:.6f} | {row.v4:.6f} | {row.v5:.6f} | {row.change_from_v4:+.6f} | {contribution} |')
    lines += ['', 'Five of six scaled components improve; direction fidelity changes by '
        f"**{record['metrics']['fid'] - float(scaled.loc[scaled.metric == 'fid', 'v4'].iloc[0]):.6f}**. "
        f"Its contribution is **{record['fidelity_overall_contribution_change_from_v4']:.6f}**, while the other five contribute "
        f"**{record['other_components_overall_contribution_change_from_v4']:+.6f}**, leaving an overall change of "
        f"**{record['overall_change_from_v4']:.6f}**. This is an arithmetic decomposition of the six-component mean, not causal attribution.", '',
        '| Raw metric | V4 | V5 | V5 minus V4 | Improvement |', '|---|---:|---:|---:|---|']
    for row in raw.itertuples():
        lines.append(f'| {row.metric} | {row.v4:.6f} | {row.v5:.6f} | {row.change_from_v4:+.6f} | {"yes" if row.improved_from_v4 else "no"} |')
    lines += ['', f"Raw expression error falls **{-record['raw_expression_error_relative_change_from_v4']:.2%}**. "
        f"Direction-fidelity yield falls **{-record['raw_fidelity_yield_relative_change_from_v4']:.2%}**. "
        'Fidelity yield is response-budget-normalized; it is not simple sign precision. '
        'The supplied aggregate metrics cannot distinguish context-specific failure, response undercoverage or incorrect directions.', '',
        'V5 simultaneously changes public H1 coverage, B/C response calibration, moment matching, control templates and A promoter corrections. '
        'No component or context can be blamed from the aggregate result. A useful follow-up is controlled source-domain ablation of these changes '
        'against V4, measuring mean/pooled effects, variance, response coverage and direction yield together. '
        'The new 300-gene neural models reconstruct control cells; their reconstruction gains do not establish better CRISPRi response predictions.', '',
        f"Rank was **1100** in the user report; the recorded API check returned **{record['rank_at_api_check']}**. "
        'Ranks depend on the current leaderboard and are not directly comparable across publication dates. '
        f"Submitted at `{record['submitted_utc']}`; reported on 2026-10-07. "
        f"Partition `{record['partition']}`, panel `{record['panel']}`, anchors `{record['anchor_set']}`.", '',
        '## Records and reproduction', '',
        'The original `submission_ready.json` remains the generation-time snapshot. Publication is recorded separately in '
        '`leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv` '
        'under `experiments/v5_atlas_calibrated/analysis/`. The local package hash is inherited from the verified readiness record; '
        'the status API does not independently expose or confirm that hash.', '',
        '```cmd', 'cd /d L:\\promec\\Animesh\\Download\\Virtual',
        '.\\.competition_env\\Scripts\\python.exe -m vcc status ' + ENTRY,
        '.\\.competition_env\\Scripts\\python.exe record_v5_result.py', '```', '',
        'Both commands read the existing entry; they do not upload a duplicate. '
        '[V5_SUBMISSION.md](V5_SUBMISSION.md) and [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md) preserve generation and verification details.', '']
    (ROOT / 'V5_LEADERBOARD_RESULTS.md').write_text('\n'.join(lines), encoding='utf-8')
    path = ROOT / 'V5_SUBMISSION.md'
    text = path.read_text(encoding='utf-8')
    text = text.replace('# V5: verified submission ready for upload', '# V5: verified and published submission')
    begin, end = '<!-- V5 publication begin -->', '<!-- V5 publication end -->'
    if begin in text:
        start = text.index(begin); stop = text.index(end, start) + len(end)
        text = text[:start] + text[stop:].lstrip('\n')
    old_intro = 'Created for the 2026 validation panel `vcc2026-val-1`. The candidate is **generated, officially packaged and independently verified**, and has not been uploaded or scored. V4 remains published with overall -0.05067391016092112; local V5 diagnostics do not establish a leaderboard improvement.'
    text = text.replace(old_intro, 'Created for the 2026 validation panel `vcc2026-val-1`. Generation, official packaging and independent verification passed. The completed leaderboard result is recorded above.')
    text = text.replace('## Submit using Windows Command Prompt', '## Original upload command (historical)')
    text = text.replace('The command uploads a new competition entry and waits for scoring. It has not been run during generation. Preserve the resulting entry ID, panel/anchor set and all six scaled/raw metrics to compare with V4. No token needs to be included in the command or committed.',
        'This is the original upload command, retained for reproduction. The entry has completed scoring; use `vcc status ' + ENTRY + '` to inspect it. Repeating the upload command creates another entry. No token belongs in a command or committed file.')
    block = begin + '\n\n' + f'Published as `{ENTRY}`, overall **{record["overall"]:.6f}**, rank 1100 at the user-reported check. V4 remains best at **-0.050674**. V5 improves five scaled components but loses direction fidelity; see [V5_LEADERBOARD_RESULTS.md](V5_LEADERBOARD_RESULTS.md) for exact comparisons and raw metrics.' + '\n\n' + end + '\n\n'
    title, rest = text.split('\n', 1)
    path.write_text(title + '\n\n' + block + rest.lstrip('\n'), encoding='utf-8')
    path = ROOT / 'README.md'; text = path.read_text(encoding='utf-8')
    text = text.replace('Three submissions (first, V2 and V4) have been published.', 'Four submissions (first, V2, V4 and V5) have been published.')
    start = text.index('## V5 atlas-calibrated submission\n'); stop = text.index('\n## ', start + 3)
    section = text[start:stop]
    section = section.replace('It has not been uploaded or scored.', f'Published as `{ENTRY}` with overall **{record["overall"]:.6f}**; V4 remains best at **-0.050674**.')
    section = section.replace('V4 remains the best published project result until V5 is scored.', 'V4 remains the best of four published project submissions. V5 improves five scaled components, but direction fidelity falls sharply. [V5_LEADERBOARD_RESULTS.md](V5_LEADERBOARD_RESULTS.md) records all scores, raw metrics and the arithmetic decomposition.')
    text = text[:start] + section + text[stop:]
    text = text.replace('V4 is the best of the three published project submissions on this validation panel.', 'V4 is the best of the four published project submissions on this validation panel; the later V5 result is documented above.')
    path.write_text(text, encoding='utf-8')
    provenance = ANALYSIS / 'documentation_manifest.json'
    if provenance.exists():
        from build_submission import sha256
        manifest = json.loads(provenance.read_text(encoding='utf-8'))
        manifest['root_code_sha256']['document_v5_submission.py'] = sha256(ROOT / 'document_v5_submission.py')
        manifest['root_code_sha256']['record_v5_result.py'] = sha256(Path(__file__))
        manifest['document_sha256']['V5_SUBMISSION.md'] = sha256(ROOT / 'V5_SUBMISSION.md')
        manifest['document_sha256']['V5_LEADERBOARD_RESULTS.md'] = sha256(ROOT / 'V5_LEADERBOARD_RESULTS.md')
        for name in ['leaderboard_result.json', 'published_api_status.json', 'leaderboard_comparison.csv', 'raw_metric_comparison.csv']:
            manifest['artifact_sha256'][name] = sha256(ANALYSIS / name)
        manifest['scope'] = 'Generation and verification records plus separately API-confirmed published V5 scores.'
        write_json(provenance, manifest)


def main():
    process = subprocess.run([sys.executable, '-m', 'vcc', '--json', 'status', ENTRY],
        capture_output=True, text=True, encoding='utf-8', check=True, timeout=45)
    api = json.loads(process.stdout)
    assert api['entry_id'] == ENTRY and api['status'] == 'published' and api['is_terminal']
    for name, value in REPORTED.items():
        np.testing.assert_allclose(api['score_avg' if name == 'overall' else 'score_' + name], value, rtol=0, atol=1e-12)
    paths = {'first': ROOT / 'analysis/leaderboard_result.json',
        'v2': ROOT / 'experiments/v2_cd4_hybrid/analysis/leaderboard_result.json',
        'v4': ROOT / 'experiments/v4_lineage_response/analysis/leaderboard_result.json'}
    previous = {key: json.loads(path.read_text(encoding='utf-8')) for key, path in paths.items()}
    for old in previous.values():
        assert api['partition'] == old['partition'] and api['panel_id'] == old['panel'] and api['anchor_version'] == old['anchor_set']
    v4 = previous['v4']; comparison = []
    for name, value in REPORTED.items():
        row = {'metric': name, **{key: old['overall'] if name == 'overall' else old['metrics'][name] for key, old in previous.items()}, 'v5': value}
        row['change_from_v4'] = value - row['v4']
        row['overall_contribution_change_from_v4'] = None if name == 'overall' else row['change_from_v4'] / 6
        comparison.append(row)
    comparison = pd.DataFrame(comparison)
    delta = REPORTED['overall'] - v4['overall']
    np.testing.assert_allclose(comparison.iloc[1:].overall_contribution_change_from_v4.sum(), delta, atol=1e-12, rtol=0)
    np.testing.assert_allclose(sum(REPORTED[k] for k in REPORTED if k != 'overall') / 6, REPORTED['overall'], atol=1e-12, rtol=0)
    raw = pd.DataFrame([{'metric': name, 'v4': v4['raw_metrics'][name], 'v5': api[name],
        'change_from_v4': api[name] - v4['raw_metrics'][name], 'higher_is_better': higher,
        'improved_from_v4': api[name] > v4['raw_metrics'][name] if higher else api[name] < v4['raw_metrics'][name]} for name, higher in RAW.items()])
    ready = json.loads((ANALYSIS / 'submission_ready.json').read_text(encoding='utf-8'))
    fid = (REPORTED['fid'] - v4['metrics']['fid']) / 6
    result = {'source': 'User-reported CLI status, independently confirmed with read-only vcc --json status.',
        'reported_date': '2026-10-07', 'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'entry_id': ENTRY, 'model_name': api['model_name'], 'final_status': 'published',
        'submitted_utc': api['submission_date'], 'rank_at_user_report': 1100, 'rank_at_api_check': api['rank'],
        'overall': REPORTED['overall'], 'metrics': {k: v for k, v in REPORTED.items() if k != 'overall'},
        'raw_metrics': {k: api[k] for k in RAW}, 'partition': api['partition'], 'panel': api['panel_id'],
        'anchor_set': api['anchor_version'], 'package_sha256': ready['sha256'],
        'package_hash_scope': 'Verified local generation record; status API does not expose a package hash.',
        'overall_change_from_v4': delta, 'fidelity_overall_contribution_change_from_v4': fid,
        'other_components_overall_contribution_change_from_v4': delta - fid,
        'raw_expression_error_relative_change_from_v4': api['expr_mse_unbiased_capped_norm'] / v4['raw_metrics']['expr_mse_unbiased_capped_norm'] - 1,
        'raw_fidelity_yield_relative_change_from_v4': api['de_wilcoxon_direction_fidelity_yield_raw'] / v4['raw_metrics']['de_wilcoxon_direction_fidelity_yield_raw'] - 1,
        'interpretation': 'Five scaled components improve; direction-fidelity decline outweighs them. V4 remains best. Aggregate metrics do not isolate causes or contexts.',
        'readiness_record_scope': 'Original submission_ready.json is preserved as the generation-time snapshot.'}
    write_json(ANALYSIS / 'published_api_status.json', api)
    write_json(ANALYSIS / 'leaderboard_result.json', result)
    comparison.to_csv(ANALYSIS / 'leaderboard_comparison.csv', index=False)
    raw.to_csv(ANALYSIS / 'raw_metric_comparison.csv', index=False)
    document_publication()
    print(comparison.to_string(index=False))
    print(raw.to_string(index=False))


if __name__ == '__main__':
    main()
