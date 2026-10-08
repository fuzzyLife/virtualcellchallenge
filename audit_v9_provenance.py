"""Check linked selection, safety and source-record hashes without modifying frozen models."""
import json
from v9_methods import ROOT,A
import build_submission as build


def main():
    selection=json.loads((A/'count_selection.json').read_text())
    confirmation=json.loads((A/'confirmation_result.json').read_text())
    model=json.loads((A/'reference_model.json').read_text())
    plan=json.loads((A/'validation_plan.json').read_text())
    neutral=json.loads((A/'decoder_neutral_manifest.json').read_text());checked=[]
    def check(path,expected):
        assert build.sha256(path)==expected,str(path)
        checked.append(str(path.relative_to(ROOT)))
    for filename,key in [('validation_plan.json','plan_sha256'),('count_calibration_scores.csv','calibration_scores_sha256'),
                         ('count_calibration_summary.csv','calibration_summary_sha256'),('decoder_neutral_manifest.json','neutral_manifest_sha256')]:
        check(A/filename,selection[key])
    for filename,key in [('select_v9_generator.py','script_sha256'),('v9_methods.py','methods_sha256'),('evaluate_v9_submission.py','evaluation_script_sha256')]:
        check(ROOT/filename,selection[key])
    check(A/'count_selection.json',confirmation['frozen_calibration_selection_sha256'])
    check(A/'count_confirmation_summary.csv',confirmation['summary_sha256'])
    for filename,key in [('decoder_neutral_scores.csv','scores_sha256'),('decoder_neutral_summary.csv','summary_sha256')]:
        check(A/filename,neutral[key])
    check(ROOT/'benchmark_v9_neutral.py',neutral['script_sha256']);check(ROOT/'v9_methods.py',neutral['methods_sha256'])
    for record in neutral['source_inputs']:check(ROOT/record['file'],record['sha256'])
    check(ROOT/'prepare_v9_validation.py',plan['script_sha256'])
    check(ROOT/'experiments/v8_corrected_transfer/analysis/validation_plan.json',plan['v8_plan_sha256'])
    check(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz',plan['source_cache_sha256'])
    check(ROOT/'experiments/v8_corrected_transfer/analysis/reference_model.json',model['prior_frozen_model_sha256'])
    check(A/'decoder_neutral_manifest.json',model['decoder_neutral_manifest_sha256'])
    for filename,digest in model['dependencies'].items():check(ROOT/filename,digest)
    used=set(plan['excluded_previous_v8_targets'])
    old=json.loads((ROOT/'experiments/v8_corrected_transfer/analysis/validation_plan.json').read_text())
    assert used=={v['target'] for c in old['contexts'] for v in c['cells']}
    for context in plan['contexts']:
        prior=json.loads((ROOT/f'experiments/v7_multisource_adaptive/analysis/cell_split_{context["family"]}.json').read_text())
        assert not {v['target'] for v in context['cells']}&(used|set(prior['targets']))
    build.write_json(A/'linked_provenance_validation.json',{'passed':True,'checked_files':sorted(set(checked)),
        'all_calibration_confirmation_neutral_source_links_match':True,'all_v7_and_v8_count_targets_excluded':True,
        'script_sha256':build.sha256(ROOT/'audit_v9_provenance.py'),
        'model_sha256':build.sha256(A/'reference_model.json'),
        'scope':'Read-only linked record integrity; no numerical model or decision changed.'})
    print('V9 linked provenance passed:',len(set(checked)),'files',flush=True)

if __name__=='__main__':main()
