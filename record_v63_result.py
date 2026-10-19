"""Record user-reported V63 publication without altering frozen experiments."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v63_identity_guard/analysis'
BASELINE = ROOT / 'experiments/v24_cd4_response/analysis/leaderboard_result.json'
PREVIOUS = ROOT / 'experiments/v61_pooled_state_response/analysis/leaderboard_result.json'
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
SCORES = dict(pds=0.6517575033476736, mse=0.09847018046403415,
              nmae=0.08629026340781953, fid=-0.01899826535007117,
              reach=0.06988104358292134, jac=0.015136044265488424)
OVERALL = 0.15042279495297764
PANEL = ('val', 'vcc2026-val-1', 'vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4')


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def close(a, b):
    assert math.isclose(a, b, rel_tol=0, abs_tol=1e-12), (a, b)


def main(check=False):
    baseline, previous = read(BASELINE), read(PREVIOUS)
    close(OVERALL, sum(SCORES.values()) / 6)
    for receipt in (baseline, previous):
        close(receipt['overall'], sum(receipt['metrics'][k] for k in METRICS) / 6)
        assert tuple(receipt[k] for k in ('partition', 'panel', 'anchor_set')) == PANEL
    rows = [dict(metric=k, v24=baseline['metrics'][k], v61=previous['metrics'][k],
                 v63=SCORES[k], change_from_v24=SCORES[k] - baseline['metrics'][k],
                 overall_contribution=(SCORES[k] - baseline['metrics'][k]) / 6)
            for k in METRICS]
    delta = OVERALL - baseline['overall']
    close(sum(r['overall_contribution'] for r in rows), delta)
    assert [r['metric'] for r in rows if r['change_from_v24'] > 0] == ['pds', 'mse', 'nmae', 'jac']
    assert [r['metric'] for r in rows if r['change_from_v24'] < 0] == ['fid', 'reach']
    ready = read(A / 'submission_ready.json')
    assert ready['experimental'] and not ready['recommended_over_V24']
    assert read(A / 'independent_validation.json')['passed']
    inputs = [BASELINE, PREVIOUS, A / 'release_manifest.json', A / 'submission_ready.json',
              A / 'independent_validation.json', ROOT / 'V63_RESULTS.md',
              ROOT / 'V63_SUBMISSION.md', Path(__file__)]
    inputs_sha = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    outputs = [A / 'leaderboard_result.json', A / 'leaderboard_comparison.csv',
               ROOT / 'V63_LEADERBOARD_RESULTS.md']
    audit_path = A / 'leaderboard_publication_audit.json'
    if check:
        audit = read(audit_path)
        assert audit['inputs_sha256'] == inputs_sha
        for name, digest in audit['outputs_sha256'].items():
            assert sha(ROOT / name) == digest, name
        result = read(outputs[0])
        assert result['metrics'] == SCORES and result['overall'] == OVERALL
        close(result['overall_change_from_v24'], delta)
        assert result['package_sha256'] == ready['sha256']
        print('V63 publication verified: same-panel arithmetic, frozen input hashes and three output hashes.')
        return
    result = dict(source='Exact published CLI scores supplied by the user; no independent API query performed.',
                  entry_id='mhKv7sajwD0DKQ8oHPYP', model_name='Identity-guarded assay response v63 experimental',
                  final_status='published', rank_at_user_report=516, overall=OVERALL, metrics=SCORES,
                  raw_metrics=None, submitted_utc=None, partition=PANEL[0], panel=PANEL[1], anchor_set=PANEL[2],
                  scoring_job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-mhkv7sajwd0dkq8ohpyp',
                  package_sha256=ready['sha256'], upload_bytes=ready['bytes'], upload_duration_seconds=984,
                  artifact_link_evidence='The user supplied the submission command naming the audited local prediction.vcc. The CLI receipt supplies no uploaded-file checksum.',
                  comparison_entry_id=baseline['entry_id'], overall_change_from_v24=delta,
                  relative_overall_change_from_v24=delta / baseline['overall'],
                  overall_change_from_v61=OVERALL - previous['overall'], recommended_over_v24=True,
                  public_development_gates_passed=False,
                  interpretation='New best published project reference on the same validation panel. Four scaled metrics improve over V24; fidelity and reach worsen. Expression is the main gain. Aggregate scores do not identify contexts or establish final-test generalization.')
    write(outputs[0], result)
    with outputs[1].open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    table = ['| Metric | V24 | V63 | Change | Contribution to overall change |',
             '| --- | ---: | ---: | ---: | ---: |']
    table += [f"| {r['metric']} | {r['v24']:.9f} | {r['v63']:.9f} | {r['change_from_v24']:+.9f} | {r['overall_contribution']:+.9f} |" for r in rows]
    report = f'''# V63 published leaderboard result

V63 published as `mhKv7sajwD0DKQ8oHPYP`, overall **{OVERALL:.9f}**, compared with V24 **{baseline['overall']:.9f}**: **{delta:+.9f} ({delta / baseline['overall']:+.2%})**. It is now the best published project reference, superseding V24 on this validation panel. User-reported rank **516** is a historical snapshot. V61 scored {previous['overall']:.9f}; V63 gains {OVERALL - previous['overall']:+.9f} over it.

## Same-panel comparison

{chr(10).join(table)}

All six values are scaled, higher-is-better accuracies, including MSE and NMAE. Their equal-weight averages reproduce the reported overalls; each change divided by six gives its contribution to the overall gain. Four metrics improve and two worsen. Raw V63 errors were not supplied: the higher expression score does not establish that raw error halved.

Expression contributes **{rows[1]['overall_contribution']:+.9f}** to overall, exceeding the net gain because fidelity and reach losses offset part of the improvements. PDS, log-fold-change accuracy and significance overlap also improve. Fidelity remains negative, and reach falls {rows[4]['change_from_v24']:.9f}; these remain clear gaps relative to V24.

## What this supports

V63 changes the source-response estimator: it removes source-common response components, translates effects through an assay-scale expression representation, and limits the correction using V24's predicted target margins. It retains V24's count generator. The published result supports retaining this complete candidate as our current validation reference. It does not isolate the benefit of any individual component; that requires component ablations with comparable validation.

The public paired profile error ratio was 0.968979699 (3.10% lower plug-in error), consistent with the direction of the hosted expression improvement. However, public count discrimination worsened while hosted PDS improved. Reach improved in development but worsened in additional public cases and in this submission. The strict public gates failed before upload and remain failed. This hosted win does not retroactively validate them: they are imperfect selection proxies, with limited unique targets and reused datasets. See the frozen [public results](V63_RESULTS.md) and [validation audit](V62_V64_VALIDATION_AUDIT.md).

## Next decision

Use V63 as the published comparison baseline, keeping V24 as the prior reference. Preserve the expression improvement while testing fidelity/reach corrections; require paired, study-held-out component ablations and shuffled controls before attributing a gain to a mechanism. Check predicted means and generated counts separately, since the public expression and discrimination checks disagree in their ability to forecast the hosted result. More complex ensembles are not justified by this result alone.

The receipt provides no per-context or per-target breakdown, so it cannot identify which context or gene caused the improvement. One aggregate validation outcome provides no estimate of uncertainty or guarantee for the final test panel.

## Provenance and reproduction

Source: the user's complete published CLI receipt, without an additional API query. Partition `{PANEL[0]}`, panel `{PANEL[1]}`, anchors `{PANEL[2]}`, matching both saved comparison receipts. Submission time was not supplied. Local package SHA-256: `{ready['sha256']}`; {ready['bytes']:,} bytes; upload duration 16 minutes 24 seconds. The user's command links this path to the entry, but the server receipt supplies no artifact hash.

The original scientific scripts, model, package, pre-upload reports and 522-file release manifest remain unchanged. Their uploaded=false and recommendation fields record pre-upload knowledge. Separate publication files record this subsequent result. README's newest status block is the current overview; older sections are historical.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe record_v63_result.py
.\\.competition_env\\Scripts\\python.exe record_v63_result.py --check
.\\.competition_env\\Scripts\\python.exe seal_v63_release.py --check
```

The recorder writes only publication files and the newest marked README block. It does not change predictions, upload files or push Git. Its check mode verifies recorded arithmetic and hashes; the release-seal check separately verifies the original scientific files.
'''
    outputs[2].write_text(report, encoding='utf-8')
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8-sig')
    begin, end = '<!-- v63-identity-guard:begin -->', '<!-- v63-identity-guard:end -->'
    assert text.count(begin) == text.count(end) == 1
    block = f'''{begin}
## V63 published: new best project score

V63 is published as `mhKv7sajwD0DKQ8oHPYP`, overall **{OVERALL:.9f}**, versus V24 **{baseline['overall']:.9f}**: **{delta:+.9f} ({delta / baseline['overall']:+.2%})**. Expression accuracy is the main gain; PDS, log-fold-change accuracy and overlap also improve. Fidelity and reach decline. **V63 is now the published reference.** Rank 516 is the user's snapshot. Public development gates remain failed; this single validation result does not establish final-test performance.

[Published scores and interpretation](V63_LEADERBOARD_RESULTS.md) · [Frozen pre-upload results](V63_RESULTS.md) · [Reproduce](V63_IDENTITY_GUARD_STRATEGY.md) · [Validation audit](V62_V64_VALIDATION_AUDIT.md). Original model, package and 522-file scientific seal are preserved. Older status blocks below describe earlier decisions.
{end}'''
    start = text.index(begin)
    finish = text.index(end, start) + len(end)
    readme.write_text(text[:start] + block + text[finish:], encoding='utf-8')
    write(audit_path, dict(source='User-reported published V63 CLI receipt', inputs_sha256=inputs_sha,
                           outputs_sha256={str(p.relative_to(ROOT)): sha(p) for p in outputs},
                           overall_reconciled=True, same_panel_and_anchors=True,
                           pre_submission_manifest_unchanged=True))
    print(f'Recorded V63 {OVERALL:.9f}; V24 delta {delta:+.9f} ({delta / baseline["overall"]:+.2%}).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify publication files without writing.')
    main(parser.parse_args().check)
