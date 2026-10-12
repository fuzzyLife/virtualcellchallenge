"""Record V16 publication while preserving frozen generation and validation evidence."""
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v16_control_preserving/analysis'
ENTRY = 'wG3GGbaTsl2dXFn7cQy1'
SCORES = dict(overall=.1416138568460195, pds=.6415529345408622,
              mse=.046085071692100965, nmae=.08034481237466624,
              fid=-.017295587840992322, reach=.08438107714423089,
              jac=.01461483316524908)
RAW = {'pds_cosine': True, 'expr_mse_unbiased_capped_norm': False,
       'de_wilcoxon_lfc_nmae': False, 'de_wilcoxon_direction_fidelity_yield_raw': True,
       'de_wilcoxon_direction_reach_raw': True, 'de_wilcoxon_sig_jaccard': True}
ARCHIVES = {'V16_RESULTS.md': 'generation_results_document.md',
            'V16_SUBMISSION.md': 'generation_submission_document.md'}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def check_close(a, b):
    assert math.isclose(a, b, rel_tol=0, abs_tol=1e-12), (a, b)


def table(rows):
    keys = list(rows[0])
    lines = ['| ' + ' | '.join(keys) + ' |', '| ' + ' | '.join(['---'] * len(keys)) + ' |']
    for row in rows:
        values = ['' if v is None else f'{v:.6f}' if isinstance(v, float) else str(v) for v in row.values()]
        lines.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join(lines) + '\n'


def save_csv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify_generation(manifest):
    for name, digest in manifest['inputs'].items():
        assert sha(ROOT / name) == digest, 'Changed generation input: ' + name
    for name, digest in manifest['outputs'].items():
        archive = A / ARCHIVES[name]
        path = archive if archive.exists() else ROOT / name
        assert sha(path) == digest, 'Changed original report: ' + name
    assert sha(ROOT / 'document_v16_submission.py') == manifest['generator_sha256']
    model = read(A / 'reference_model.json')
    for name, digest in model['dependencies'].items():
        assert sha(ROOT / name) == digest, 'Changed model dependency: ' + name
    assert (model['candidate'], model['response'], model['emitter']) == ('residual10', 'v13', 'residual10')
    return len(model['dependencies'])


def main(api_file=None):
    manifest = read(A / 'report_manifest.json')
    seal_hash = sha(A / 'report_manifest.json')
    dependency_count = verify_generation(manifest)
    if api_file:
        api = read(api_file)
    else:
        query = subprocess.run([sys.executable, '-m', 'vcc', '--json', 'status', ENTRY],
                               capture_output=True, text=True, encoding='utf-8', check=True, timeout=45)
        api = json.loads(query.stdout)
    assert api['entry_id'] == ENTRY and api['status'] == 'published' and api['is_terminal']
    for metric, value in SCORES.items():
        check_close(api['score_avg' if metric == 'overall' else 'score_' + metric], value)
    baseline_path = ROOT / 'experiments/v13_response_programs/analysis/leaderboard_result.json'
    baseline = read(baseline_path)
    assert (api['partition'], api['panel_id'], api['anchor_version']) == (baseline['partition'], baseline['panel'], baseline['anchor_set'])
    delta = SCORES['overall'] - baseline['overall']
    scaled = []
    for metric, value in SCORES.items():
        before = baseline['overall'] if metric == 'overall' else baseline['metrics'][metric]
        scaled.append(dict(metric=metric, v13=before, v16=value, change_from_v13=value-before,
                           overall_contribution_change=None if metric == 'overall' else (value-before)/6))
    check_close(sum(v for k, v in SCORES.items() if k != 'overall')/6, SCORES['overall'])
    check_close(sum(row['overall_contribution_change'] or 0 for row in scaled), delta)
    assert sum(SCORES[k] > baseline['metrics'][k] for k in RAW_SCORE_KEYS) == 4
    raw = [dict(metric=k, v13=baseline['raw_metrics'][k], v16=api[k],
                change_from_v13=api[k]-baseline['raw_metrics'][k],
                relative_change=api[k]/baseline['raw_metrics'][k]-1,
                higher_is_better=higher,
                improved=api[k] > baseline['raw_metrics'][k] if higher else api[k] < baseline['raw_metrics'][k])
           for k, higher in RAW.items()]
    ready, audit = read(A/'submission_ready.json'), read(A/'independent_validation.json')
    assert ready['official_prep_passed'] and audit['passed']
    assert ready['sha256'] == audit['package_sha256']
    assert (A.parent/'prediction.vcc').stat().st_size == ready['bytes']
    previous = [read(p)['overall'] for p in ROOT.glob('experiments/v*/analysis/leaderboard_result.json') if p.parent != A]
    previous.append(read(ROOT/'analysis/leaderboard_result.json')['overall'])
    best_so_far = SCORES['overall'] > max(previous)
    assert best_so_far, 'Best-score wording requires V16 to exceed every other recorded submission.'
    result = dict(source='User-reported CLI publication confirmed by read-only vcc --json status.',
                  recorded_utc=datetime.now(timezone.utc).isoformat(), entry_id=ENTRY,
                  model_name=api['model_name'], final_status='published', submitted_utc=api['submission_date'],
                  rank_at_user_report=562, rank_at_api_check=api['rank'], overall=SCORES['overall'],
                  metrics={k:v for k,v in SCORES.items() if k!='overall'},
                  raw_metrics={k:api[k] for k in RAW}, partition=api['partition'], panel=api['panel_id'],
                  anchor_set=api['anchor_version'], package_sha256=ready['sha256'], upload_bytes=ready['bytes'],
                  upload_duration_seconds=986, overall_change_from_v13=delta,
                  relative_overall_change_from_v13=delta/baseline['overall'],
                  report_manifest_sha256=seal_hash, recorder_sha256=sha(__file__),
                  best_recorded_project_validation_score=best_so_far,
                  interpretation='Four scaled components improve; fidelity and NMAE decline. V13 response arrays are retained; count generation changes. A single published result does not establish significance or final-round generalization.')
    write_json(A/'published_api_status.json', api)
    write_json(A/'leaderboard_result.json', result)
    save_csv(A/'leaderboard_comparison.csv', scaled)
    save_csv(A/'raw_metric_comparison.csv', raw)
    for name, archived in ARCHIVES.items():
        path = A/archived
        if not path.exists():
            assert sha(ROOT/name) == manifest['outputs'][name]
            path.write_bytes((ROOT/name).read_bytes())
        assert sha(path) == manifest['outputs'][name]
    block = (f"V16 is published as `{ENTRY}`, overall **{SCORES['overall']:.9f}**, a gain of **{delta:+.9f}** "
             f"({delta/baseline['overall']:.2%}) over V13 on the same panel and anchors. It is the best recorded project validation score. "
             f"Rank was **562** in the user report and **{api['rank']}** at the API check; historical ranks are not directly comparable.")
    report = f'''# V16 published leaderboard result

{block}

## Scaled comparison

{table(scaled)}
All six displayed components are scaled scores: higher is better, including MSE and NMAE. Four improve; fidelity and fold-change accuracy decline. Expression accuracy contributes +0.002201 and reach +0.001468 to overall. Contributions are score changes divided by six, reconcile with the overall gain, and are arithmetic rather than causal attribution.

## Raw comparison

{table(raw)}
Raw expression error decreases {-(api['expr_mse_unbiased_capped_norm']/baseline['raw_metrics']['expr_mse_unbiased_capped_norm']-1):.2%}; raw reach increases {(api['de_wilcoxon_direction_reach_raw']/baseline['raw_metrics']['de_wilcoxon_direction_reach_raw']-1):.2%}. Fidelity remains below its zero-score anchor. Improved count generation has not resolved response direction or biological dispersion.

## What we learned

Production retains the exact V13 response arrays and changes the population generator to residual10: 0.25 expected + 0.65 Poisson(expected) + 0.10 empirical response counts, followed by the original dual-moment solver. Gene-calibrated shrinkage and balanced molecule transport failed the frozen development guards and were not deployed. This result supports retaining residual10 as the new reference, with only modest evidence from one leaderboard experiment.

Local reserved comparisons improved raw NMAE by 0.00304047 on the eligible H1 panels, but the hosted NMAE score declines. Missing HCT116/HEK293T NMAE panels prevented complete six-metric verification. Fidelity also declines on the reserved and hosted comparisons. Local tests remain exploratory on partly reused biology; the later score does not retroactively turn them into independent or complete validation.

Residual10 retains only 20.5% of real total variance in the public control diagnostic and still produces substantial false DE there. Preserving control variance with balanced transport fixed that diagnostic but worsened actual perturbation accuracy. Next tests should hold residual10 fixed while separating response-direction errors from distribution errors with real treatment-replica controls, eligible per-family metrics, and matched seed/mean ablations. Another broad ensemble is not supported by these results: the V15 ensemble and V16 shrinkage alternatives were rejected.

## Provenance and reproduction

Submitted `{api['submission_date']}`. Upload took 16 minutes 26 seconds. Partition `{api['partition']}`, panel `{api['panel_id']}`, anchors `{api['anchor_version']}`. Package: 4,073,615,360 bytes; SHA-256 `{ready['sha256']}`. Official prep and independent full-payload/provenance/replay validation passed before upload. The API does not expose the uploaded file hash; the user's upload record connects the audited local artifact to this entry. This recorder checks file size and the prior sealed validation records, without rehashing the 3.79 GiB artifact.

The original report manifest, its {len(manifest['inputs'])} inputs, all {dependency_count} model dependencies, generation-time uploaded=false snapshots, and numerical code remain unchanged. Both original report outputs are archived byte-for-byte. `publication_manifest.json` maps their original hashes to those archives and separately records the current publication documents. After rerunning the generation-time report writer, rerun this recorder to restore the publication overlay. README is a project index and is not included in the immutable report hashes.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
.\\.competition_env\\Scripts\\python.exe record_v16_result.py --api-file experiments/v16_control_preserving/analysis/published_api_status.json
```
'''
    (ROOT/'V16_LEADERBOARD_RESULTS.md').write_text(report, encoding='utf-8', newline='\n')
    original = (A/ARCHIVES['V16_RESULTS.md']).read_text(encoding='utf-8')
    original = original.replace('It has not been uploaded; no score or rank gain is known.', 'At generation time it had not been uploaded and no score was known; publication is recorded below.')
    original = original.replace('V13 remains the best published project score (0.13788313368012328).', 'V13 was the best published project score at generation time (0.13788313368012328).')
    original = original.replace('No upload is performed.', 'The generation pipeline performs no upload; the user subsequently uploaded this package.')
    title, body = original.split('\n', 1)
    (ROOT/'V16_RESULTS.md').write_text(title+'\n\n'+block+' See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md).\n\n'+body.lstrip('\n'), encoding='utf-8', newline='\n')
    submission = f'''# V16 published submission

{block} See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md) for all components, raw comparisons and limitations.

- File: `experiments\\v16_control_preserving\\prediction.vcc`.
- Candidate: `residual10`; response: exact V13; generator: residual10.
- Size: 4,073,615,360 bytes (3.79 GiB).
- SHA-256: `{ready['sha256']}`.
- Shape: 360,000 cells x 18,533 genes; 900 groups, exactly 400 cells each.
- Official prep and independent complete-payload/provenance/replay validation passed.
- Available public guards passed; complete six-metric public verification was unavailable.

Inspect the already published entry from Command Prompt:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {ENTRY}
```

The completed upload used model name `V13 response residual dispersion v16`. Repeating submit would create a duplicate experiment. Original generation-time upload instructions are archived in `experiments/v16_control_preserving/analysis/generation_submission_document.md`. [V16_RESULTS.md](V16_RESULTS.md) preserves the actual public tests and their limitations; the publication overlay does not alter their outcomes.
'''
    (ROOT/'V16_SUBMISSION.md').write_text(submission, encoding='utf-8', newline='\n')
    readme = ROOT/'README.md'
    text = readme.read_text(encoding='utf-8')
    text = text.replace('Ten submissions (first, V2, V4, V5, V6, V7, V9, V10, V11 and V13) have been published.', 'Eleven submissions (first, V2, V4, V5, V6, V7, V9, V10, V11, V13 and V16) have been published.')
    text = text.replace('V13 is the best recorded validation submission, with overall score **0.137883134**, improving on V11 by **0.000879160** (0.64%).', f"V16 is the best recorded validation submission, with overall score **{SCORES['overall']:.9f}**, improving on V13 by **{delta:.9f}** ({delta/baseline['overall']:.2%}).")
    text = text.replace('See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md) for the current result and [V13_SUBMISSION.md](V13_SUBMISSION.md) for reproducibility.', 'See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md) for the current result and [V16_SUBMISSION.md](V16_SUBMISSION.md) for reproducibility.')
    text = text.replace('V13 is now the best published reference.', 'V16 is now the best published reference.').replace('V13 is now best.', 'V16 is now best.')
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith('The latest and best published submission is'):
            lines[i] = f"The latest and best published submission is V16 (`{ENTRY}`), overall **{SCORES['overall']:.9f}**. See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md)."
        if line.startswith('Published as `AvAsYkuxFSjbxYDCd0wJ`'):
            lines[i] = line.replace('This is the best recorded project validation score.', 'This was the best recorded project validation score at V13 publication; V16 now improves it.')
    text = '\n'.join(lines)+'\n'
    start, end = '<!-- v16-research:begin -->', '<!-- v16-research:end -->'
    assert start in text and end in text
    before, tail = text.split(start, 1)
    _, after = tail.split(end, 1)
    section = f'''\n\n## V16 published result

{block} Four scaled components improve, led by expression accuracy and reach; fidelity and NMAE decline. The exact V13 response predictor is retained, with residual10 count generation.

Frozen development compares 252 populations (224 new, 28 verified reused baseline), plus 18 control diagnostics and 36 reserved comparisons. Calibrated response shrinkage and balanced transport were rejected. Available public guards passed, with missing NMAE panels limiting complete validation. Official packaging and independent full-payload/provenance/replay checks passed. Methods: [V16_STRATEGY.md](V16_STRATEGY.md). Public tests: [V16_RESULTS.md](V16_RESULTS.md). Hosted comparison: [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md). Artifact/status: [V16_SUBMISSION.md](V16_SUBMISSION.md).\n\n'''
    readme.write_text(before+start+section+end+after, encoding='utf-8', newline='\n')
    verify_generation(manifest)
    assert sha(A/'report_manifest.json') == seal_hash
    artifacts = ['published_api_status.json', 'leaderboard_result.json', 'leaderboard_comparison.csv',
                 'raw_metric_comparison.csv', *ARCHIVES.values()]
    write_json(A/'publication_manifest.json', dict(entry_id=ENTRY, recorder_sha256=sha(__file__),
               publication_artifact_sha256={n:sha(A/n) for n in artifacts},
               publication_document_sha256={n:sha(ROOT/n) for n in ['V16_RESULTS.md', 'V16_SUBMISSION.md', 'V16_LEADERBOARD_RESULTS.md']},
               report_manifest_sha256=seal_hash,
               original_report_output_overlay={n:str((A/p).relative_to(ROOT)).replace('\\','/') for n,p in ARCHIVES.items()},
               comparison_baseline_sha256=sha(baseline_path),
               all_generation_file_hashes_verified_before_and_after=True,
               verified_model_dependencies=dependency_count,
               scope='Publication overlay only; frozen predictor, local acceptance, generation evidence and package unchanged.'))
    print(table(scaled))
    print('V16 publication recorded; frozen generation evidence unchanged.')


RAW_SCORE_KEYS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--api-file', type=Path)
    main(parser.parse_args().api_file)
