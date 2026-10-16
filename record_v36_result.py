"""Record user-reported published V36 scores without changing sealed experiment evidence."""
from pathlib import Path
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
A = ROOT / 'experiments/v36_population_views/analysis'
SCORES = dict(overall=0.14184548169072003, pds=0.6433442472382521,
              mse=0.0501735761598077, nmae=0.08093425364417435,
              fid=-0.01643177767448657, reach=0.07816022080205248,
              jac=0.014892369974520037)

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    previous = read(ROOT / 'experiments/v24_cd4_response/analysis/leaderboard_result.json')
    ready = read(A / 'submission_ready.json')
    assert ready['official_prep_passed']
    assert (A.parent / 'prediction.vcc').stat().st_size == ready['bytes']
    assert read(A / 'final_delivery_audit.json')['passed']
    assert math.isclose(SCORES['overall'], sum(SCORES[k] for k in SCORES if k != 'overall') / 6, abs_tol=1e-12)
    panel = ('val', 'vcc2026-val-1', 'vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4')
    assert panel == (previous['partition'], previous['panel'], previous['anchor_set'])
    before = dict(overall=previous['overall'], **previous['metrics'])
    changes = {k: SCORES[k] - before[k] for k in SCORES}
    assert math.isclose(changes['overall'], sum(changes[k] for k in changes if k != 'overall') / 6, abs_tol=1e-12)
    sealed = {p: sha(p) for p in A.rglob('*.json') if p.name != 'leaderboard_result.json'}
    result = {
        'source': 'User-provided published vcc CLI output; no independent API confirmation or remote file checksum claimed.',
        'entry_id': 'GftPza37V9QL6u7omFOD',
        'model_name': 'Retrieval-contrast response correction v36',
        'final_status': 'published', 'rank_at_user_report': 543,
        'overall': SCORES['overall'],
        'metrics': {k: v for k, v in SCORES.items() if k != 'overall'},
        'partition': panel[0], 'panel': panel[1], 'anchor_set': panel[2],
        'package_sha256': ready['sha256'], 'upload_bytes': ready['bytes'],
        'upload_duration_seconds': 973,
        'package_link_provenance': 'User upload command links entry to local package; checksum copied from audited generation receipt.',
        'generation_ready_sha256': sha(A / 'submission_ready.json'),
        'scaled_scores_higher_is_better': True,
        'comparison_entry_id': previous['entry_id'],
        'changes_from_v24': changes,
        'overall_contributions_vs_v24': {k: v / 6 for k, v in changes.items() if k != 'overall'},
        'best_recorded_project_validation_score': False,
        'preferred_submission': 'experiments/v24_cd4_response/prediction.vcc',
        'interpretation': 'Small regression against V24 on the same panel and anchors. Expression accuracy and DE overlap improve; PDS, NMAE, FID and especially REACH decline. This tests the selected quarter-strength retrieval contrast, not the rejected full population-view model. One aggregate result does not establish statistical significance or final-test performance.'
    }
    path = A / 'leaderboard_result.json'
    if path.exists():
        assert read(path) == result, 'Existing published receipt differs; refusing overwrite.'
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')
    import document_v35_v36_results
    document_v35_v36_results.main()
    assert all(sha(p) == digest for p, digest in sealed.items()), 'Sealed experiment evidence changed.'
    print(json.dumps({'entry': result['entry_id'], 'overall_change': changes['overall'], 'sealed_json_unchanged': len(sealed)}, indent=2))

if __name__ == '__main__':
    main()
