"""Record V69's supplied published receipt; preserve sealed scientific artifacts."""
import argparse
import csv
from decimal import Decimal as D
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/v69_published'
SCIENCE = ROOT / 'experiments/v69_count_bridge/analysis'
BASE = ROOT / 'experiments/v66_assay_strength/analysis/leaderboard_result.json'
REPORT = ROOT / 'V69_LEADERBOARD_RESULTS.md'
VALUES = dict(overall='0.09526284481907588', pds='0.4722661732603027',
              mse='0.00966231669808748', nmae='0.08325305971139878',
              jac='0.008059922552236526')
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
BEGIN, END = '<!-- v69-published:begin -->', '<!-- v69-published:end -->'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    path.write_bytes((json.dumps(value, indent=2, allow_nan=False) + '\n').encode('utf-8'))


def main(check=False):
    baseline = read(BASE)
    ready = read(SCIENCE / 'submission_ready.json')
    validation = read(SCIENCE / 'independent_validation.json')
    confirmation = read(SCIENCE / 'confirmation_assessment.json')
    manifest = read(SCIENCE / 'release_manifest.json')
    assert validation['passed'] and validation['package_sha256'] == ready['sha256']
    assert ready['experimental'] and not ready['recommended_over_V66']
    assert not confirmation['release_supported']
    panel = dict(partition='val', panel='vcc2026-val-1',
                 anchor_set='vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4')
    assert all(baseline[k] == v for k, v in panel.items())
    assert abs(sum(D(str(baseline['metrics'][m])) for m in METRICS) / 6
               - D(str(baseline['overall']))) < D('1e-15')
    overall = D(VALUES['overall'])
    delta = overall - D(str(baseline['overall']))
    relative = delta / D(str(baseline['overall']))
    rows = []
    for metric in METRICS:
        old = D(str(baseline['metrics'][metric]))
        new = D(VALUES[metric]) if metric in VALUES else None
        rows.append(dict(metric=metric, v66=str(old), v69=str(new) if new is not None else None,
                         change_from_v66=str(new-old) if new is not None else None,
                         status='reported' if new is not None else 'not reported'))
    known_delta = sum(D(r['change_from_v66']) for r in rows if r['v69'] is not None)
    ratios = confirmation['ratios']['hybrid25']
    error_ratio = sum(ratios.values()) / len(ratios)
    inputs = [BASE, Path(__file__), *(SCIENCE / p for p in (
        'submission_ready.json', 'independent_validation.json',
        'confirmation_assessment.json', 'release_manifest.json'))]
    hashes = {p.relative_to(ROOT).as_posix(): sha(p) for p in inputs}
    paths = [OUT / 'leaderboard_result.json', OUT / 'leaderboard_comparison.csv', REPORT]
    audit_path = OUT / 'publication_audit.json'
    result = dict(source='Published CLI receipt supplied by user; no independent API query.',
                  entry_id='dj3mrGx3Md59leld1ZI3',
                  model_name='V66 probability-coupled count mixture v69 experimental',
                  final_status='published', rank_at_user_report=772,
                  overall=float(overall), exact_reported_values=VALUES,
                  metrics={m: float(VALUES[m]) if m in VALUES else None for m in METRICS},
                  unreported_metrics=['fid', 'reach'], raw_metrics=None,
                  submitted_utc=None, **panel,
                  scoring_job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-dj3mrgx3md59leld1zi3',
                  local_package_sha256=ready['sha256'], local_package_bytes=ready['bytes'],
                  artifact_link_evidence='User command identifies audited local package; server receipt has no checksum.',
                  upload_duration_seconds=910, comparison_entry_id=baseline['entry_id'],
                  overall_change_from_v66=str(delta), relative_change_from_v66=str(relative),
                  best_published_version='V66', public_release_supported=False,
                  public_mean_study_profile_error_ratio=error_ratio)
    table = '\n'.join(f"| {r['metric'].upper()} | {float(r['v66']):.6f} | "
                      + (f"{float(r['v69']):.6f} | {float(r['change_from_v66']):+.6f} |"
                         if r['v69'] is not None else 'unreported | unavailable |') for r in rows)
    report = f'''# V69 published result: count mixture rejected

Entry `{result['entry_id']}` is **published**, according to the user-supplied CLI receipt. Overall **{overall:.9f}** versus V66 **{baseline['overall']:.9f}**: **{delta:+.9f} ({relative:+.2%})**. V66 remains the highest published project score. Rank 772 is a snapshot; ranks from different dates are not a controlled comparison.

| Scaled score (higher is better) | V66 | V69 | Change |
|---|---:|---:|---:|
| Overall | {baseline['overall']:.6f} | {overall:.6f} | {delta:+.6f} |
{table}

FID and REACH are absent from the supplied receipt, so their individual values remain null. MSE and NMAE above are scaled scores, not raw errors. The four reported component declines contribute **{known_delta / 6:+.9f}** to overall under the established six-metric equal-weight rule. PDS and expression accuracy account for most of the drop.

## What the experiment establishes

V69 changed 100 of 400 cell positions with a probability-coupled count bridge while retaining 300 exact V66 rows. Frozen V66 response means were generator inputs; this did not guarantee preservation of realized expression profiles. Official preparation and independent full-payload validity checks passed. A valid submission can still predict the biology poorly.

Public confirmation had already rejected the candidate: mean paired study profile-error ratio **{error_ratio:.6f}** (about **{(error_ratio-1)*100:.2f}% worse**), with failed fidelity, reach and Jaccard guards and unavailable aggregate NMAE. Public PDS increased from **{confirmation['summary']['v66']['pds']:.6f}** to **{confirmation['summary']['hybrid25']['pds']:.6f}**, but hosted PDS declined. Beating the matched shuffled mixture was insufficient evidence of a useful change. These reused public studies do not reliably rank candidates for hidden contexts.

This rejects this particular quarter-population bridge as an improvement over V66. It does not isolate one cause or disprove all distribution models. Mean calibration, sampling, source-assay mismatch and metric transfer are confounded. Another blend or a larger variance change is not supported by this result.

## Next decision

Retain V66 as the reference. Before another release, require improvement over V66 on held-out response profiles and discrimination, with adequate DE measurements and representative control nulls. Prioritize an assay-aware response-mean/count-likelihood teacher on an independent study or reserved donor; compare against V66 and shuffled evidence before generating cells. Freeze the split and thresholds before opening that reserve. This is a proposed experiment, not a completed model or a promised ranking gain.

## Reproduction and provenance

Run `python record_submission_result_69.py`, then `python record_submission_result_69.py --check` from the root. Only the standard library is required. JSON and CSV comparisons are in `analysis/v69_published/`; input/output hashes are in `publication_audit.json`. Exact reported decimal strings are retained.

The original scientific release remains sealed: {len(manifest['files'])} files, local package SHA256 `{ready['sha256']}`. The user's command links the uploaded artifact; the receipt supplies no server checksum. Historical pre-upload files still say not uploaded and remain unchanged; this publication record supersedes that status. The recorder performs no API query, training, prediction generation, upload or GitHub push.
'''
    block = f'''{BEGIN}
## V69 published: count mixture regressed; V66 remains best

V69 entry `{result['entry_id']}` scored **{overall:.9f}**, versus V66 **{baseline['overall']:.9f}**: **{delta:+.9f} ({relative:+.2%})**. All four reported component scores declined, led by discrimination and expression accuracy. FID and REACH were not included in the supplied receipt. Public confirmation guards had failed, and its PDS gain reversed on the leaderboard. **Retain V66 as the best published reference.**

[Published scores and lessons](V69_LEADERBOARD_RESULTS.md) · [Reproducible result recorder](record_submission_result_69.py). Original sealed artifacts are preserved. Older blocks below describe historical pre-upload decisions.
{END}'''
    readme_path = ROOT / 'README.md'
    readme = readme_path.read_bytes().decode('utf-8-sig')
    if check:
        audit = read(audit_path)
        assert audit['inputs_sha256'] == hashes
        for name, digest in audit['outputs_sha256'].items():
            assert sha(ROOT / name) == digest, name
        assert read(paths[0]) == result
        assert REPORT.read_text(encoding='utf-8') == report
        assert readme.count(BEGIN) == readme.count(END) == 1 and block in readme
        with paths[1].open(encoding='utf-8', newline='') as stream:
            assert list(csv.DictReader(stream)) == [{k: '' if v is None else v for k, v in r.items()} for r in rows]
        print('V69 publication verified: exact scores, missing metrics, arithmetic and source/output hashes.')
        return
    OUT.mkdir(parents=True, exist_ok=True)
    write(paths[0], result)
    with paths[1].open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    REPORT.write_bytes(report.encode('utf-8'))
    assert readme.count(BEGIN) == readme.count(END)
    if BEGIN in readme:
        assert readme.count(BEGIN) == 1
        start, finish = readme.index(BEGIN), readme.index(END) + len(END)
        readme = readme[:start] + block + readme[finish:]
    else:
        readme = block + '\n\n' + readme
    readme_path.write_bytes(readme.encode('utf-8'))
    write(audit_path, dict(inputs_sha256=hashes,
                          outputs_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
                          same_panel_and_anchors=True, unreported_metrics_preserved_as_null=True))
    print(f'Recorded V69 {overall:.9f}; V66 change {delta:+.9f} ({relative:+.2%}).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    main(parser.parse_args().check)
