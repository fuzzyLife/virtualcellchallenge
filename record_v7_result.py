"""Record V7 publication using read-only status and preserve generation snapshots."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import build_submission as build
from record_v5_result import RAW

ROOT = build.ROOT
A = ROOT / 'experiments/v7_multisource_adaptive/analysis'
ENTRY = 'k0pkYuJRcw3C0U8K0NgL'
REPORTED = {'overall': -0.05138675064550524, 'pds': 0.23585266664934945,
            'mse': 0., 'nmae': 0.07186092628425385, 'fid': -0.6726758462891342,
            'reach': 0.0867146928111978, 'jac': -0.030072943328698278}


def overlay(text, block):
    begin, end = '<!-- V7 publication begin -->', '<!-- V7 publication end -->'
    if begin in text:
        start = text.index(begin)
        stop = text.index(end, start) + len(end)
        text = text[:start] + text[stop:].lstrip('\n')
    title, rest = text.split('\n', 1)
    return title + '\n\n' + begin + '\n\n' + block + '\n\n' + end + '\n\n' + rest.lstrip('\n')


def document_publication():
    result = json.loads((A / 'leaderboard_result.json').read_text(encoding='utf-8'))
    scaled = pd.read_csv(A / 'leaderboard_comparison.csv')
    raw = pd.read_csv(A / 'raw_metric_comparison.csv')
    lines = ['# V7 published leaderboard result', '',
        f'Entry `{ENTRY}` is **published**, overall **{result["overall"]:.9f}**, rank {result["rank_at_user_report"]} in the user report. '
        'V4 remains best at **-0.050673910**. All six published project entries share the same validation panel and anchors.', '',
        '| Scaled metric | V4 | V6 | V7 | V7 minus V4 | Overall contribution |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for row in scaled.itertuples():
        contribution = '' if row.metric == 'overall' else f'{row.overall_contribution_change_from_v4:+.9f}'
        lines.append(f'| {row.metric} | {row.v4:.9f} | {row.v6:.9f} | {row.v7:.9f} | {row.change_from_v4:+.9f} | {contribution} |')
    lines += ['', f'Overall is **{result["overall_change_from_v4"]:+.9f}** versus V4 and **{result["overall_change_from_v6"]:+.9f}** versus V6. '
        'Fidelity and LFC accuracy improve; discrimination, reach and overlap decline; scaled expression remains zero. '
        'The discrimination loss offsets the fidelity and LFC gains. These contributions are arithmetic attribution within the six-component mean.', '',
        '| Raw metric | V4 | V6 | V7 | V7 minus V4 | Improved versus V4 |',
        '| --- | ---: | ---: | ---: | ---: | --- |']
    for row in raw.itertuples():
        lines.append(f'| {row.metric} | {row.v4:.9f} | {row.v6:.9f} | {row.v7:.9f} | {row.change_from_v4:+.9f} | {"yes" if row.improved_from_v4 else "no"} |')
    lines += ['', f'Raw expression error changes **{result["raw_expression_error_relative_change_from_v4"]:+.3%}**, '
        f'while raw fidelity yield changes **{result["raw_fidelity_yield_relative_change_from_v4"]:+.3%}**. '
        'The zero scaled expression score is a clipped score, not zero prediction error.', '',
        '## What this experiment establishes', '',
        'The V7 and V4 response-effect arrays match exactly. V7 selected the V4 transfer fallback because none of the expanded-source, pooling, neural or ridge alternatives passed the fixed calibration guards. '
        'Count generation uses the same zero-aware binomial/Poisson recipe with a new seed and different guide/depth-stratified control draws. '
        'This publication is a resampling comparison for the existing estimator; it does not measure the hidden performance of the rejected alternatives.', '',
        'The overall difference is small. Two realizations cannot estimate a reliable seed-variation distribution or establish statistical significance. '
        'Aggregate scores do not identify the weak context or target. Keep V4 as the published reference and diagnose source-to-context response calibration and unsupported DE calls before adopting a more complex model. '
        'Require improvements across paired offline seeds; another unchanged-model upload is not a biological modeling advance.', '',
        '## Provenance and reproduction', '',
        f'Submitted `{result["submitted_utc"]}`; upload took 13 minutes 18 seconds. Rank {result["rank_at_api_check"]} at the recorded API check; ranks can move. '
        f'Partition `{result["partition"]}`, panel `{result["panel"]}`, anchors `{result["anchor_set"]}`.', '',
        'Generation-time readiness, model, package, independent-validation and pipeline records are preserved. Publication is stored separately in `leaderboard_result.json`, '
        '`published_api_status.json`, `leaderboard_comparison.csv`, `raw_metric_comparison.csv` and `publication_manifest.json`. '
        'The package SHA-256 identifies the independently checked local artifact; the status API does not expose the uploaded package hash.', '',
        '```bat', r'cd /d L:\promec\Animesh\Download\Virtual',
        r'.\.competition_env\Scripts\python.exe -m vcc status ' + ENTRY,
        r'.\.competition_env\Scripts\python.exe record_v7_result.py', '```', '',
        'These inspect the existing entry. [V7_SUBMISSION.md](V7_SUBMISSION.md) preserves the original model selection, validation and historical upload command.', '']
    (ROOT / 'V7_LEADERBOARD_RESULTS.md').write_text('\n'.join(lines), encoding='utf-8')
    path = ROOT / 'V7_SUBMISSION.md'
    text = path.read_text(encoding='utf-8')
    block = (f'Published as `{ENTRY}`, overall **{result["overall"]:.9f}**, rank {result["rank_at_user_report"]} in the user report. '
             'V4 remains best. This is the same response estimator with new control sampling/count draws. '
             'See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md) for API-confirmed scaled/raw comparisons.')
    text = overlay(text, block)
    text = text.replace('The upload-ready file is', 'The verified package is')
    text = text.replace('The package has not been uploaded and has no hidden leaderboard score.', 'The package is published; its leaderboard result is recorded above.')
    text = text.replace('## Upload using Windows CMD', '## Historical upload command (already published)')
    text = text.replace('If waiting is interrupted after upload, resume status polling for the returned entry ID rather than submitting again:',
        'The command above created the published entry. Repeating it creates a duplicate submission. Inspect the existing result with:')
    text = text.replace('vcc status ENTRY_ID --wait', 'vcc status ' + ENTRY)
    path.write_text(text, encoding='utf-8')
    path = ROOT / 'README.md'
    text = path.read_text(encoding='utf-8')
    text = text.replace('Five submissions (first, V2, V4, V5 and V6) have been published.', 'Six submissions (first, V2, V4, V5, V6 and V7) have been published.')
    rows = text.splitlines()
    for index, row in enumerate(rows):
        if row.startswith('The latest published submission is'):
            rows[index] = f'The latest published submission is V7 (`{ENTRY}`), with overall **{result["overall"]:.9f}**. V4 remains best. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md) and [V7_SUBMISSION.md](V7_SUBMISSION.md).'
        if row.startswith(('V7 is ready for upload at', 'V7 is published as')):
            rows[index] = (f'V7 is published as `{ENTRY}`, overall **{result["overall"]:.9f}** versus V4 **-0.050673910**. '
                'It uses exactly V4 response effects with new count/control sampling. Expanded HCT116/HEK293T sources and neural/ridge alternatives were tested, but none passed the calibration guards. '
                'Official packaging and the independent full-payload audit passed. V4 remains best. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md), '
                '[V7_SUBMISSION.md](V7_SUBMISSION.md) and [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md).')
    text = '\n'.join(rows) + '\n'
    text = text.replace('best of five published project submissions', 'best of six published project submissions').replace('best of the five published project submissions', 'best of the six published project submissions')
    path.write_text(text, encoding='utf-8')
    path = ROOT / 'V7_MODEL_STRATEGY.md'
    rows = path.read_text(encoding='utf-8').splitlines()
    rows[2] = f'Status: independently verified and published as `{ENTRY}`, overall **{result["overall"]:.9f}**. The response arrays equal V4 exactly; the new count realization did not surpass V4. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md) for the publication and [V7_SUBMISSION.md](V7_SUBMISSION.md) for selection, hashes and reproduction.'
    path.write_text('\n'.join(rows) + '\n', encoding='utf-8')
    path = ROOT / 'PERFORMANCE_IMPROVEMENT_PLAN.md'
    text = path.read_text(encoding='utf-8').replace('It has not been uploaded.', f'It is published as `{ENTRY}`, overall **{result["overall"]:.9f}**, slightly below V4. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md).')
    path.write_text(text, encoding='utf-8')
    path = ROOT / 'analysis/v7_source_inventory/manifest_scope.json'
    scope = json.loads(path.read_text(encoding='utf-8'))
    scope.update(status='Completed, independently verified and published; V4 remains best', uploaded=True,
                 publication={'entry_id': ENTRY, 'overall': result['overall'], 'source': 'Read-only API-confirmed publication record'})
    build.write_json(path, scope)
    names = ['leaderboard_result.json', 'published_api_status.json', 'leaderboard_comparison.csv', 'raw_metric_comparison.csv']
    build.write_json(A / 'publication_manifest.json', {
        'entry_id': ENTRY, 'recorder_sha256': build.sha256(Path(__file__)),
        'report_generator_sha256': build.sha256(ROOT / 'document_v7_submission.py'),
        'publication_artifact_sha256': {name: build.sha256(A / name) for name in names},
        'publication_document_sha256': build.sha256(ROOT / 'V7_LEADERBOARD_RESULTS.md'),
        'submission_document_sha256': build.sha256(ROOT / 'V7_SUBMISSION.md'),
        'generation_snapshot_sha256': {name: build.sha256(A / name) for name in ['submission_ready.json', 'independent_validation.json', 'reference_model.json', 'pipeline_status.json']},
        'scope': 'Publication overlay. Generation snapshots retain their historical upload=false state and original document hashes.'})


def main(api_file=None):
    if api_file:
        data = Path(api_file).read_bytes()
        api = json.loads(data.decode('utf-16' if data.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig'))
    else:
        process = subprocess.run([sys.executable, '-m', 'vcc', '--json', 'status', ENTRY],
            capture_output=True, text=True, encoding='utf-8', check=True, timeout=45)
        api = json.loads(process.stdout)
    assert api['entry_id'] == ENTRY and api['status'] == 'published' and api['is_terminal']
    for key, value in REPORTED.items():
        np.testing.assert_allclose(api['score_avg' if key == 'overall' else 'score_' + key], value, rtol=0, atol=1e-12)
    previous = {name: json.loads((ROOT / directory / 'leaderboard_result.json').read_text(encoding='utf-8'))
        for name, directory in {'first': Path('analysis'), 'v2': Path('experiments/v2_cd4_hybrid/analysis'),
        'v4': Path('experiments/v4_lineage_response/analysis'), 'v5': Path('experiments/v5_atlas_calibrated/analysis'),
        'v6': Path('experiments/v6_fidelity_calibrated/analysis')}.items()}
    for old in previous.values():
        assert (api['partition'], api['panel_id'], api['anchor_version']) == (old['partition'], old['panel'], old['anchor_set'])
    scaled = []
    for key, value in REPORTED.items():
        row = {'metric': key, **{name: old['overall'] if key == 'overall' else old['metrics'][key] for name, old in previous.items()}, 'v7': value}
        row.update(change_from_v4=value-row['v4'], change_from_v6=value-row['v6'],
                   overall_contribution_change_from_v4=None if key == 'overall' else (value-row['v4'])/6)
        scaled.append(row)
    scaled = pd.DataFrame(scaled)
    delta = REPORTED['overall'] - previous['v4']['overall']
    np.testing.assert_allclose(scaled.iloc[1:].overall_contribution_change_from_v4.sum(), delta, rtol=0, atol=1e-12)
    np.testing.assert_allclose(np.mean([v for k, v in REPORTED.items() if k != 'overall']), REPORTED['overall'], rtol=0, atol=1e-12)
    raw = pd.DataFrame([{'metric': key, 'v4': previous['v4']['raw_metrics'][key], 'v6': previous['v6']['raw_metrics'][key], 'v7': api[key],
        'change_from_v4': api[key]-previous['v4']['raw_metrics'][key], 'higher_is_better': higher,
        'improved_from_v4': api[key] > previous['v4']['raw_metrics'][key] if higher else api[key] < previous['v4']['raw_metrics'][key]} for key, higher in RAW.items()])
    ready = json.loads((A / 'submission_ready.json').read_text(encoding='utf-8'))
    validation = json.loads((A / 'independent_validation.json').read_text(encoding='utf-8'))
    equality = json.loads((A / 'v4_fallback_effect_comparison.json').read_text(encoding='utf-8'))
    assert validation['passed'] and validation['package_sha256'] == ready['sha256']
    assert equality['maximum_absolute_log2_difference'] == 0 and equality['coordinates_different_above_1e_6'] == 0
    assert equality['v7_model_arrays_sha256'] == build.sha256(A / 'reference_model_arrays.npz')
    result = {'source': 'User-reported CLI publication, independently confirmed by read-only vcc --json status.',
        'recorded_utc': datetime.now(timezone.utc).isoformat(), 'entry_id': ENTRY, 'model_name': api['model_name'],
        'final_status': 'published', 'submitted_utc': api['submission_date'], 'rank_at_user_report': 1052, 'rank_at_api_check': api['rank'],
        'overall': REPORTED['overall'], 'metrics': {k: v for k, v in REPORTED.items() if k != 'overall'}, 'raw_metrics': {k: api[k] for k in RAW},
        'partition': api['partition'], 'panel': api['panel_id'], 'anchor_set': api['anchor_version'],
        'package_sha256': ready['sha256'], 'upload_bytes': ready['bytes'], 'upload_duration_seconds': 798,
        'package_hash_scope': 'Independently validated local artifact; status API does not expose the uploaded package hash.',
        'overall_change_from_v4': delta, 'overall_change_from_v6': REPORTED['overall']-previous['v6']['overall'],
        'raw_expression_error_relative_change_from_v4': api['expr_mse_unbiased_capped_norm']/previous['v4']['raw_metrics']['expr_mse_unbiased_capped_norm']-1,
        'raw_fidelity_yield_relative_change_from_v4': api['de_wilcoxon_direction_fidelity_yield_raw']/previous['v4']['raw_metrics']['de_wilcoxon_direction_fidelity_yield_raw']-1,
        'response_effect_arrays_equal_v4': True, 'interpretation': 'V4 remains best. Mean effects are unchanged; new control draws/count sampling give a slightly lower overall result. This does not evaluate the rejected new response estimators or establish a seed uncertainty distribution.',
        'readiness_record_scope': 'Generation-time snapshots are preserved; this separate publication record gives current scoring evidence.',
        'recorder_sha256': build.sha256(Path(__file__))}
    build.write_json(A / 'published_api_status.json', api)
    build.write_json(A / 'leaderboard_result.json', result)
    scaled.to_csv(A / 'leaderboard_comparison.csv', index=False)
    raw.to_csv(A / 'raw_metric_comparison.csv', index=False)
    document_publication()
    print(scaled[['metric', 'v4', 'v7', 'change_from_v4']].to_string(index=False))
    print('API-confirmed V7 publication documented; generation artifacts unchanged.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--api-file', type=Path, help='Captured read-only vcc --json status output')
    main(parser.parse_args().api_file)
