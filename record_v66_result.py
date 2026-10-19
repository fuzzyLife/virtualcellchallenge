"""Record the user's published V66 receipt without changing sealed experiments."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v66_assay_strength/analysis'
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
SCORES = dict(pds=0.6517009514165119, mse=0.09657486167945828,
              nmae=0.08527061361861961, fid=-0.02260859548432163,
              reach=0.07962563498187135, jac=0.014729659086930079)
OVERALL = 0.15088218754984492
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
    p63 = ROOT / 'experiments/v63_identity_guard/analysis/leaderboard_result.json'
    p24 = ROOT / 'experiments/v24_cd4_response/analysis/leaderboard_result.json'
    v63, v24 = read(p63), read(p24)
    close(OVERALL, sum(SCORES.values()) / 6)
    for baseline in (v63, v24):
        assert tuple(baseline[k] for k in ('partition', 'panel', 'anchor_set')) == PANEL
        close(baseline['overall'], sum(baseline['metrics'][k] for k in METRICS) / 6)
    local = read(A / 'count_assessment.json')
    sensitivity = read(A / 'paired_profile_sensitivity.json')
    ready = read(A / 'submission_ready.json')
    assert ready['experimental'] and not ready['recommended_over_V63']
    assert read(A / 'independent_validation.json')['passed']
    assert not local['release_supported'] and not local['checks']['reach']
    mapping = dict(pds=('pds', 1), mse=('profile_square_error', -1),
                   nmae=('nmae', -1), fid=('fidelity', 1),
                   reach=('reach', 1), jac=('jaccard', 1))
    rows = []
    for metric in METRICS:
        field, orientation = mapping[metric]
        local_delta = (local['summary']['assay75'][field] - local['summary']['v63'][field]) * orientation
        hosted_delta = SCORES[metric] - v63['metrics'][metric]
        rows.append(dict(metric=metric, v24=v24['metrics'][metric],
                         v63=v63['metrics'][metric], v66=SCORES[metric],
                         change_from_v63=hosted_delta,
                         overall_contribution_from_v63=hosted_delta / 6,
                         change_from_v24=SCORES[metric] - v24['metrics'][metric],
                         local_field=field, local_change_oriented_higher_better=local_delta,
                         local_hosted_direction_agrees=(local_delta * hosted_delta > 0)))
    delta = OVERALL - v63['overall']
    close(sum(r['overall_contribution_from_v63'] for r in rows), delta)
    assert [r['metric'] for r in rows if r['change_from_v63'] > 0] == ['reach']
    assert [r['metric'] for r in rows if r['local_hosted_direction_agrees']] == ['jac']
    inputs = [p63, p24, A / 'release_manifest.json', A / 'submission_ready.json',
              A / 'independent_validation.json', A / 'count_assessment.json',
              A / 'paired_profile_sensitivity.json',
              ROOT / 'experiments/v65_closed_response/analysis/broader_finalist_diagnostic.json',
              ROOT / 'V66_RESULTS.md', ROOT / 'V65_V66_SUBMISSION.md', Path(__file__)]
    inputs_sha = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    outputs = [A / 'leaderboard_result.json', A / 'leaderboard_comparison.csv',
               ROOT / 'V66_LEADERBOARD_RESULTS.md']
    audit_path = A / 'leaderboard_publication_audit.json'
    if check:
        audit = read(audit_path)
        assert audit['inputs_sha256'] == inputs_sha
        for name, digest in audit['outputs_sha256'].items():
            assert sha(ROOT / name) == digest, name
        result = read(outputs[0])
        assert result['metrics'] == SCORES and result['overall'] == OVERALL
        assert tuple(result[k] for k in ('partition', 'panel', 'anchor_set')) == PANEL
        assert result['entry_id'] == 'x8V3fUSgFIoAzHbNx4eX'
        assert result['package_sha256'] == ready['sha256']
        close(result['overall_change_from_v63'], delta)
        assert result['local_hosted_direction_disagreements'] == 5
        print('V66 publication verified: same-panel arithmetic, five direction disagreements and all input/output hashes.')
        return
    result = dict(source='Exact published CLI receipt supplied by user; no independent API query.',
                  entry_id='x8V3fUSgFIoAzHbNx4eX', model_name='Assay-strength response v66 experimental',
                  final_status='published', rank_at_user_report=518, overall=OVERALL, metrics=SCORES,
                  raw_metrics=None, submitted_utc=None, partition=PANEL[0], panel=PANEL[1], anchor_set=PANEL[2],
                  scoring_job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-x8v3fusgfioazhbnx4ex',
                  package_sha256=ready['sha256'], upload_bytes=ready['bytes'], upload_duration_seconds=978,
                  artifact_link_evidence='User command names the audited local V66 package; server receipt supplies no artifact checksum.',
                  comparison_entry_id=v63['entry_id'], overall_change_from_v63=delta,
                  relative_overall_change_from_v63=delta / v63['overall'],
                  overall_change_from_v24=OVERALL - v24['overall'],
                  highest_published_project_score=True, statistically_established_improvement=False,
                  public_development_gates_passed=False, local_hosted_direction_disagreements=5,
                  interpretation='Tiny new highest validation score, entirely reach-driven; five other scaled metrics decline. Keep V63 as near-tied reference. Local count-test directions disagree for five metrics; no final-test guarantee or component-level causal attribution.')
    write(outputs[0], result)
    with outputs[1].open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    table = ['| Scaled metric | V63 | V66 | Change | Overall contribution |',
             '| --- | ---: | ---: | ---: | ---: |']
    table += [f"| {r['metric']} | {r['v63']:.9f} | {r['v66']:.9f} | {r['change_from_v63']:+.9f} | {r['overall_contribution_from_v63']:+.9f} |" for r in rows]
    directions = ['| Metric | Small public count test | Hosted validation |', '| --- | --- | --- |']
    directions += [f"| {r['metric']} | {'improves' if r['local_change_oriented_higher_better'] > 0 else 'worsens'} | {'improves' if r['change_from_v63'] > 0 else 'worsens'} |" for r in rows]
    reach_contribution = next(r['overall_contribution_from_v63'] for r in rows if r['metric'] == 'reach')
    report = f'''# V66 published result: small reach-driven gain

V66 published as `x8V3fUSgFIoAzHbNx4eX`, overall **{OVERALL:.9f}**, versus V63 **{v63['overall']:.9f}**: **{delta:+.9f} ({delta / v63['overall']:+.2%})**. This is the highest published project score, but the gain is small and its statistical significance is unknown. Keep V63 as a near-tied comparison reference.

The user-reported rank is **518**, compared with V63's earlier 516. These are different-time leaderboard snapshots: other teams' results change, so a higher score can accompany a worse rank. Rank movement alone does not measure this model change.

## Component arithmetic

{chr(10).join(table)}

All six CLI values are scaled, higher-is-better accuracies, including MSE and NMAE. Their arithmetic mean reproduces each overall. Reach is the only improvement. It contributes **{reach_contribution:+.9f}** to overall; the other five metrics together contribute **{delta - reach_contribution:+.9f}**. Their sum is the small net gain. Fidelity is more negative, and expression, discrimination, fold-change accuracy and significance overlap all decline.

Against V24, overall gains **{OVERALL - v24['overall']:+.9f} ({(OVERALL - v24['overall']) / v24['overall']:+.2%})**: PDS, expression and fold-change accuracy improve; fidelity, reach and overlap worsen. The prior 18-publication review remains a frozen historical snapshot; V66 is the nineteenth publication now recorded.

## What changed and what the tests actually predicted

V66 raises the original assay-response proposal from half strength to three-quarter strength before applying the same V24 target-margin guard. The representation, source evidence and residual count generator are inherited from V63. The guard can reduce the accepted step, so the final effects are not uniformly multiplied by 1.5. See the sealed [method and results](V66_RESULTS.md). This whole-version comparison does not isolate a causal benefit from strength alone.

{chr(10).join(directions)}

Five of six directions disagree; overlap alone agrees. Public error quantities are oriented so that smaller profile error and NMAE error count as improvements. Their magnitudes cannot be compared directly with hosted scaled accuracies. The public count test had 30 unique target/context pairs and two technical seeds; its 240 comparison records include shared baselines rather than 240 independent biological cases.

The paired profile-error ratio was **{sensitivity['ratio']:.9f}** (about 0.37% lower error), concentrated in Orion; H1 and CD4 were worse. Leave-one-target/context-out sensitivity crossed a ratio of one. The predeclared strict REACH guard failed, yet hosted REACH is the sole gain. These observations show that passing most small-cohort proxies did not forecast the hosted component changes.

The broader 655 target/context-pair desired-mean diagnostic already warned of slightly worse profile error (equal-study ratio 1.000903825) and discrimination (PDS change -0.001063193). Those directions agree with the hosted expression and PDS declines. This broader check measures intended bulk means before cell generation, and it reuses public biology; it is not an independent official-score validation set.

## Next experimental priority

Treat V66 as a narrow validation gain, not a broadly better response model. Before another strength grid or ensemble, improve the selection benchmark: compare small and broad cohorts on exactly the same target/context pairs, replay intended means and generated counts separately, and test whether study-held-out component rankings transfer to a genuinely untouched public dataset. Freeze selection rules before observing that holdout. Retain shuffled-response and unchanged-generator controls. This result motivates those tests; they have not been run as part of recording this receipt.

No per-context/per-target scores, raw errors, confidence intervals or server-side package checksum were provided. Aggregate scores cannot identify which context or gene produced the gain, quantify significance, or establish performance on the final test panel.

## Provenance and reproduction

Source: the user's complete published CLI receipt; no additional API request. Partition `{PANEL[0]}`, panel `{PANEL[1]}`, anchors `{PANEL[2]}`, identical to V63 and V24. Submission timestamp is unknown. Upload duration was 16 minutes 18 seconds. Local package: {ready['bytes']:,} bytes, SHA-256 `{ready['sha256']}`. The command links the local path to the entry; the server receipt supplies no checksum to independently confirm uploaded bytes.

The original 653-file V65/V66 scientific seal and prior 522-file V63 seal remain unchanged. Their pre-upload status/recommendation fields are historical. Separate publication records and the newest README block describe this later result. The publication audit hashes its source receipts, sealed diagnostics, recorder and three output files; README remains mutable.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe record_v66_result.py
.\\.competition_env\\Scripts\\python.exe record_v66_result.py --check
.\\.competition_env\\Scripts\\python.exe seal_v65_v66_release.py --check
```

Recording and checking do not change predictions, upload files or push Git.
'''
    outputs[2].write_text(report, encoding='utf-8')
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8-sig')
    begin, end = '<!-- v65-v66-release:begin -->', '<!-- v65-v66-release:end -->'
    assert text.count(begin) == text.count(end) == 1
    block = f'''{begin}
## V66 published: small new best score, driven by reach

V66 is published as `x8V3fUSgFIoAzHbNx4eX`, overall **{OVERALL:.9f}**, versus V63 **{v63['overall']:.9f}**: **{delta:+.9f} ({delta / v63['overall']:+.2%})**. Reach alone improves; the other five scaled metrics decline. **V66 has the highest published project score; V63 remains a near-tied reference.** Rank 518 is a time-specific snapshot. Statistical significance and final-test improvement are unknown.

The small public count test predicted the wrong direction for five of six metrics; the broader 655-pair desired-mean check had warned of expression and discrimination regressions. The strict public REACH guard remains failed. Improve validation transfer before more tuning. The 480 V65/V66 count-comparison records include audited shared baselines and repeated biology.

[Published scores, comparison and validation mismatch](V66_LEADERBOARD_RESULTS.md) · [Frozen 18-submission history](V65_SUBMISSION_HISTORY.md) · [Pre-upload V66 results](V66_RESULTS.md) · [Submission provenance](V65_V66_SUBMISSION.md). The original model, package and scientific seals are preserved. Older blocks below describe earlier decisions.
{end}'''
    start = text.index(begin)
    finish = text.index(end, start) + len(end)
    readme.write_text(text[:start] + block + text[finish:], encoding='utf-8')
    write(audit_path, dict(source='User-reported published V66 CLI receipt',
                           inputs_sha256=inputs_sha,
                           outputs_sha256={str(p.relative_to(ROOT)): sha(p) for p in outputs},
                           overall_reconciled=True, same_panel_and_anchors=True,
                           local_hosted_direction_disagreements=5,
                           original_scientific_seals_unchanged=True))
    print(f'Recorded V66 {OVERALL:.9f}; V63 delta {delta:+.9f} ({delta / v63["overall"]:+.2%}); reach alone improves.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify publication hashes and arithmetic without writing.')
    main(parser.parse_args().check)
