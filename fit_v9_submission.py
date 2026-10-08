"""Freeze unchanged V4 responses with V9's calibration-selected count decoder."""
import json,shutil
import build_submission as build
from v9_methods import ROOT,A

if __name__=='__main__':
    if (A/'reference_model.json').exists():raise FileExistsError('Production model already frozen')
    prior=ROOT/'experiments/v8_corrected_transfer/analysis'
    original=json.loads((prior/'reference_model.json').read_text())
    assert original['response']=='v4_a100'
    assert original['model_arrays_sha256']==build.sha256(prior/'reference_model_arrays.npz')
    for file,digest in original['dependencies'].items():assert build.sha256(ROOT/file)==digest
    selection=json.loads((A/'count_selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text())
    assert confirmation['frozen_calibration_selection_sha256']==build.sha256(A/'count_selection.json')
    shutil.copyfile(prior/'reference_model_arrays.npz',A/'reference_model_arrays.npz')
    files=['v9_methods.py','prepare_v9_validation.py','evaluate_v9_submission.py','select_v9_generator.py',
        'confirm_v9_selection.py','benchmark_v9_neutral.py','fit_v9_submission.py','build_v9_submission.py','validate_v9_submission.py']
    dependencies={**original['dependencies'],**{f:build.sha256(ROOT/f) for f in files}}
    record={**original,'model':'V9 zero-aware activation-preserving partial moment calibration',
        'candidate':selection['candidate'],'generator':selection['generator'],
        'confirmation_improvement':confirmation['confirmed_improvement'],'confirmation_utility_difference':confirmation['utility_difference'],
        'validation_plan_sha256':build.sha256(A/'validation_plan.json'),
        'count_selection_sha256':build.sha256(A/'count_selection.json'),
        'confirmation_result_sha256':build.sha256(A/'confirmation_result.json'),
        'decoder_neutral_manifest_sha256':build.sha256(A/'decoder_neutral_manifest.json'),
        'dependencies':dependencies,'prior_frozen_model_sha256':build.sha256(prior/'reference_model.json'),
        'source_policy':'Exact V4/V7 effects preserved; isolate initialization and moment-correction strength. H1 excluded from its own source benchmark; HCT/HEK absent from this V4 response.',
        'scope':'New public count-validation targets after V8 exploration; no hidden challenge labels. Failed confirmation is reported without reselection.'}
    record.pop('profile_selection_sha256',None)
    build.write_json(A/'reference_model.json',record)
    print('V9 production frozen',selection['candidate'],flush=True)
