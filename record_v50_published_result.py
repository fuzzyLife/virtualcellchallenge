"""Record user-reported V50 publication separately from immutable pre-upload proofs."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v50_response_prior/analysis'
KEYS = ['pds', 'mse', 'nmae', 'fid', 'reach', 'jac']
BASELINE = dict(overall=0.14299785521745675, pds=0.644246366710111,
                mse=0.04756840736075477, nmae=0.08227057284911048,
                fid=-0.01579183745269364, reach=0.08483984010194835,
                jac=0.01485378173550954)
SCORES = dict(overall=0.14152262164295007, pds=0.6445079188293643,
              mse=0.04429553851868299, nmae=0.08339644503073829,
              fid=-0.015672608176707365, reach=0.07760600943183436,
              jac=0.015002426223787849)

def read_result():
    p = A / 'published_result.json'
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None

def amend_report(text):
    r = read_result()
    if not r:
        return text
    text = text.replace('All new results are raw public diagnostics, not private anchor-scaled leaderboard scores.',
        'Experiment tables below contain raw public diagnostics. The subsequent published-result table contains official scaled scores reported by the user.')
    text = text.replace('Private leaderboard score: **unknown**.',
        'Private leaderboard score was unknown at the pre-upload decision; the published result is recorded below.')
    text = text.replace('Official packaging passed; not uploaded.',
        'Official packaging passed before upload. The user subsequently uploaded the package; see the published result below.')
    lines = ['\n## Published V50 result\n',
        'Source: user-provided VCC CLI output; this record does not claim a fresh API query. '
        'Entry `u0bZUCY7giU6jkxbuiys` is **published**, rank **557 at reporting**, overall **0.14152262164295007**. '
        'Panel `vcc2026-val-1`, partition `val`; anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.\n',
        '| Scaled score (higher is better) | V24 | V50 | V50 minus V24 |',
        '|---|---:|---:|---:|']
    for k in ['overall'] + KEYS:
        lines.append(f'| {k} | {BASELINE[k]:.9f} | {SCORES[k]:.9f} | {r["comparison"]["delta_vs_v24"][k]:+.9f} |')
    lines += ['', '**V24 remains the best published submission.** REACH and expression accuracy losses outweigh the four small gains. '
        'Historical ranks are not directly comparable because the leaderboard changes over time. All six entries above are scaled accuracy scores, including MSE and NMAE; higher is better.', '',
        'The public REACH gain failed to transfer in direction. The seed-sensitive public surrogate therefore did not reliably rank this small change. '
        'This is evidence against promoting the tested amplitude calibration over V24, not proof that every response prior or neural architecture is ineffective. '
        'Before another release, investigate this public/private ranking mismatch and require improvements across generation seeds and held-out biological datasets; D3 remains unopened.', '',
        'The frozen pre-upload receipts retain their original unknown-score and uploaded-false fields as historical snapshots. '
        'The separate [published result](experiments/v50_response_prior/analysis/published_result.json) binds this amendment to the validated package SHA256. '
        'Recreate the amendment with `SCIENCE_PYTHON record_v50_published_result.py`, then refresh documentation with `SCIENCE_PYTHON finalize_docs_v50.py`.', '']
    return text + '\n'.join(lines)

def amend_readme(text):
    if not read_result():
        return text
    text = text.replace('The surrogate gain is seed-sensitive; private superiority is unproven.',
        'The surrogate gain was seed-sensitive; the subsequently published V50 result does not beat V24.')
    text = text.replace('V24 remains best published at0.14299785521745675. New private scores are unknown;',
        'V50 is published at **0.141522622** (entry `u0bZUCY7giU6jkxbuiys`, rank 557 at reporting), **-0.001475234** versus V24. '
        'REACH and expression accuracy decline despite four small gains. V24 remains best published at **0.14299785521745675**. '
        'The public REACH improvement did not transfer; the public selection surrogate needs stronger validation. '
        'The published amendment is separate from frozen pre-upload receipts;')
    return text

def main():
    # Check reported scaled-score aggregation and the exact package identity in the ready receipt.
    for values in [BASELINE, SCORES]:
        assert abs(sum(values[k] for k in KEYS) / 6 - values['overall']) < 1e-12
    proof_names = ['submission_ready.json', 'final_release_manifest.json', 'release_decision.json']
    proofs = {n: hashlib.sha256((A / n).read_bytes()).hexdigest() for n in proof_names}
    ready = json.loads((A / 'submission_ready.json').read_text())
    assert ready['sha256'] == '4629cfca53b031527823e18792d51902fde9bd6e19e84513c4f2dcfc7945a4cb'
    result = dict(provenance='User-provided CLI output; not independently queried from API',
        entry='u0bZUCY7giU6jkxbuiys', model='Empirical response-prior calibration v50',
        status='published', rank_at_reporting=557, publication_timestamp=None,
        job='projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-u0bzucy7giu6jkxbuiys',
        partition='val', panel='vcc2026-val-1',
        anchor_set='vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4',
        upload=dict(reported_size='3.8 GiB', duration='16m14s', average_rate='4.0 MiB/s'),
        package_sha256=ready['sha256'], package_bytes=ready['bytes'],
        scaled_scores=SCORES, all_scaled_metrics_higher_is_better=True,
        comparison=dict(baseline='V24', baseline_entry='81Y2yX7utdAau9rr3jCj',
            baseline_scaled_scores=BASELINE,
            delta_vs_v24={k:SCORES[k]-BASELINE[k] for k in ['overall']+KEYS},
            champion='V24', significant_improvement_established=False),
        frozen_pre_upload_receipt_sha256=proofs)
    old = read_result()
    if old is not None:
        assert old == result, 'Existing published record differs; review rather than overwrite.'
    else:
        (A / 'published_result.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    assert proofs == {n:hashlib.sha256((A / n).read_bytes()).hexdigest() for n in proof_names}
    print(json.dumps(result['comparison'], indent=2))

if __name__ == '__main__':
    main()
