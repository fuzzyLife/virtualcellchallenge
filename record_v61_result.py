"""Record user-reported V61 publication; preserve all pre-upload scientific seals."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v61_pooled_state_response/analysis'
BASELINE = ROOT / 'experiments/v24_cd4_response/analysis/leaderboard_result.json'
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
SCORES = dict(pds=0.5562245219230239, mse=0.0,
              nmae=0.08266919718674648, fid=-0.01890699218246981,
              reach=0.06994726866634961, jac=0.010646848339915903)
OVERALL = 0.11676347398892765


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
    baseline = read(BASELINE)
    close(OVERALL, sum(SCORES.values()) / 6)
    close(baseline['overall'], sum(baseline['metrics'][k] for k in METRICS) / 6)
    panel = ('val', 'vcc2026-val-1', 'vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4')
    assert tuple(baseline[k] for k in ('partition', 'panel', 'anchor_set')) == panel
    rows = [dict(metric=k, v24=baseline['metrics'][k], v61=SCORES[k],
                 change=SCORES[k] - baseline['metrics'][k],
                 overall_contribution=(SCORES[k] - baseline['metrics'][k]) / 6)
            for k in METRICS]
    delta = OVERALL - baseline['overall']
    close(sum(r['overall_contribution'] for r in rows), delta)
    assert [r['metric'] for r in rows if r['change'] > 0] == ['nmae']
    major_share = sum(r['overall_contribution'] for r in rows if r['metric'] in ('pds', 'mse')) / delta
    inputs = [BASELINE, A / 'release_manifest.json', A / 'submission_ready.json',
              A / 'independent_validation.json', ROOT / 'V61_RESULTS.md',
              ROOT / 'V61_SUBMISSION.md', Path(__file__)]
    inputs_sha = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    outputs = [A / 'leaderboard_result.json', A / 'leaderboard_comparison.csv',
               ROOT / 'V61_LEADERBOARD_RESULTS.md']
    audit_path = A / 'leaderboard_publication_audit.json'
    if check:
        audit = read(audit_path)
        assert audit['inputs_sha256'] == inputs_sha
        for name, digest in audit['outputs_sha256'].items():
            assert sha(ROOT / name) == digest, name
        result = read(outputs[0])
        assert result['metrics'] == SCORES and result['overall'] == OVERALL
        print('V61 publication verified: same-panel arithmetic, input seals and three output hashes.')
        return
    ready = read(A / 'submission_ready.json')
    assert ready['experimental'] and not ready['recommended_over_V24']
    assert read(A / 'independent_validation.json')['passed']
    result = dict(source='Exact published CLI scores supplied by the user; no independent API query performed.',
                  entry_id='JNIsuxo0KEVQveGuZoLX', model_name='Pooled state-response transfer v61 experimental',
                  final_status='published', rank_at_user_report=674, overall=OVERALL, metrics=SCORES,
                  raw_metrics=None, submitted_utc=None, partition=panel[0], panel=panel[1], anchor_set=panel[2],
                  scoring_job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-jnisuxo0kevqveguzolx',
                  package_sha256=ready['sha256'], upload_bytes=ready['bytes'], upload_duration_seconds=987,
                  artifact_link_evidence='The user supplied the submission command naming the audited local prediction.vcc. The CLI receipt supplies no uploaded-file checksum.',
                  comparison_entry_id=baseline['entry_id'], overall_change_from_v24=delta,
                  relative_overall_change_from_v24=delta / baseline['overall'],
                  pds_and_expression_share_of_net_decline=major_share, recommended_over_v24=False,
                  interpretation='Five scaled metrics regress; only NMAE improves. The result agrees directionally with failed public transfer gates; aggregate scores cannot identify affected contexts or prove a causal mechanism.')
    write(outputs[0], result)
    with outputs[1].open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    table = ['| Metric | V24 | V61 | Change | Contribution to overall change |',
             '| --- | ---: | ---: | ---: | ---: |']
    table += [f"| {r['metric']} | {r['v24']:.9f} | {r['v61']:.9f} | {r['change']:+.9f} | {r['overall_contribution']:+.9f} |" for r in rows]
    report = f'''# V61 published leaderboard result

V61 published as `JNIsuxo0KEVQveGuZoLX`, overall **{OVERALL:.9f}**, compared with V24 **{baseline['overall']:.9f}**: **{delta:+.9f} ({delta / baseline['overall']:+.2%})**. V24 remains the recommended published reference. Rank **674** is the user's snapshot; historical ranks change as other teams submit.

## Same-panel comparison

{chr(10).join(table)}

All six components above are scaled, higher-is-better scores, including MSE and NMAE. Their equal-weight average reproduces each reported overall; component changes divided by six reproduce the overall decline. The V61 expression score of zero means it did not earn a positive scaled expression improvement over the baseline. It does not mean zero raw prediction error. Raw V61 metrics were not supplied and are not inferred.

PDS and expression together account for **{major_share:.1%}** of the net decline. NMAE is the sole gain, only {rows[2]['change']:+.9f}; fidelity, reach and overlap also worsen. This arithmetic locates the lost score, but does not isolate a causal mechanism or identify context A, B or C as responsible.

## Interpretation and next decision

The state-response approach improved within-dataset reconstruction for CD4 and H1 but failed its held-study response and count-generation support gates before upload. Its equal-study response MSE ratio was 1.025339 versus V24 (lower is better), and count-profile expression error worsened in every held study. The hosted regression agrees with that warning. See [frozen public results](V61_RESULTS.md) and [V59-V61 review](V59_V61_STATE_REVIEW.md).

The evidence does not support stronger weighting of this state correction or an ensemble that assumes it adds independent useful signal. Retain V24 as the reference. Future candidates should pass the predeclared held-study response and integer-count gates before another submission is recommended. These reused public datasets cannot establish final-test generalization, and this one aggregate leaderboard result does not prove the exact cause of failure.

## Provenance and reproduction

Source: the user's complete published CLI receipt, without an additional API query. Partition `{panel[0]}`, panel `{panel[1]}`, anchors `{panel[2]}`, matching the saved V24 receipt. No submission date was supplied. The recorded local package hash is `{ready['sha256']}` ({ready['bytes']:,} bytes); the user's command links the entry to this path, but the receipt supplies no server-side artifact hash. Upload duration: 16 minutes 27 seconds.

The original scientific code, model, package, pre-upload reports and 297-file release manifest are preserved. Their historical uploaded=false fields describe generation time. The separate publication receipt, comparison CSV and publication audit record this subsequent outcome. README is the current status overview.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe record_v61_result.py
.\\.competition_env\\Scripts\\python.exe record_v61_result.py --check
.\\.competition_env\\Scripts\\python.exe seal_v61_release.py --check
```

Reproduction writes only publication outputs and the marked README status block. It does not fit models, alter predictions, submit a package or push Git.
'''
    outputs[2].write_text(report, encoding='utf-8')
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8-sig')
    begin, end = '<!-- v61-state-release:begin -->', '<!-- v61-state-release:end -->'
    assert text.count(begin) == text.count(end) == 1
    block = f'''{begin}
## V61 published: state-response transfer regressed

V61 is published as `JNIsuxo0KEVQveGuZoLX`, overall **{OVERALL:.9f}**, versus V24 **{baseline['overall']:.9f}** ({delta / baseline['overall']:+.2%}). PDS and expression explain {major_share:.1%} of the net decline; only NMAE improves slightly. This agrees with the failed public support gates. **V24 remains recommended.** User-reported rank 674 is a historical snapshot.

[Published scores and interpretation](V61_LEADERBOARD_RESULTS.md) · [Preserved pre-upload results](V61_RESULTS.md) · [V59-V61 review](V59_V61_STATE_REVIEW.md) · [Reproduce experiments](V61_POOLED_STATE_STRATEGY.md). Original model, package, scientific code and release seal remain unchanged. No automatic upload.
{end}'''
    start = text.index(begin)
    finish = text.index(end, start) + len(end)
    readme.write_text(text[:start] + block + text[finish:], encoding='utf-8')
    write(audit_path, dict(source='User-reported published V61 CLI receipt', inputs_sha256=inputs_sha,
                           outputs_sha256={str(p.relative_to(ROOT)): sha(p) for p in outputs},
                           overall_reconciled=True, same_panel_and_anchors=True,
                           pre_submission_manifest_unchanged=True))
    print(f'Recorded V61 {OVERALL:.9f}; V24 delta {delta:+.9f} ({delta / baseline["overall"]:+.2%}).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify publication files without writing.')
    main(parser.parse_args().check)
