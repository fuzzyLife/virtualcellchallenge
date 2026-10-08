"""Verify final V10 diagnostic, production, packaging and frozen-code record links."""
import json,sys,platform
from importlib.metadata import version
from v10_methods import ROOT,A
import build_submission as build


def main():
    checked=[]
    def read(name):return json.loads((A/name).read_text(encoding='utf-8'))
    def check(path,expected):
        assert build.sha256(path)==expected,str(path)
        checked.append(str(path.relative_to(ROOT)))
    ready=read('submission_ready.json');audit=read('independent_validation.json')
    depth=read('production_depth_audit.json');corrected=read('corrected_depth_validation.json')
    model=read('reference_model.json');stage=read('memory_safe_stage.json')
    assert ready['official_prep_passed'] and audit['passed'] and depth['passed']
    assert audit['verified_rows']==depth['verified_cells']==360000
    assert audit['package_sha256']==depth['package_sha256']==ready['sha256']
    assert stage['source_sha256']==audit['prediction_sha256']
    assert (ROOT/'experiments/v10_population_transport/prediction.vcc').stat().st_size==ready['bytes']
    check(A/'reference_model.json',ready['model_sha256'])
    check(A/'reference_model.json',audit['model_sha256'])
    check(ROOT/'build_v10_submission.py',ready['generator_sha256'])
    check(ROOT/'v10_methods.py',ready['methods_sha256'])
    check(ROOT/'package_v10_memory_safe.py',ready['packaging_adapter_sha256'])
    check(ROOT/'package_v10_memory_safe.py',stage['adapter_sha256'])
    check(A/'memory_safe_stage.json',ready['memory_safe_stage_sha256'])
    check(ROOT/'validate_v10_submission.py',audit['validator_sha256'])
    check(A/'prediction_control_rows.csv',depth['provenance_sha256'])
    check(ROOT/'audit_v10_validation_depths.py',depth['script_sha256'])
    for name,key in [('count_selection.json','selection_sha256'),('validation_plan.json','plan_sha256'),
                     ('corrected_depth_scores.csv','scores_sha256'),('corrected_depth_summary.csv','summary_sha256')]:
        check(A/name,corrected[key])
    check(ROOT/'audit_v10_validation_depths.py',corrected['script_sha256'])
    check(ROOT/'v10_methods.py',corrected['methods_sha256'])
    for name,digest in model['dependencies'].items():check(ROOT/name,digest)
    uncertainty=read('paired_uncertainty.json')
    check(ROOT/'analyze_v10_tradeoffs.py',uncertainty['script_sha256'])
    check(A/'count_calibration_scores.csv',uncertainty['source_score_hashes']['calibration'])
    check(A/'corrected_depth_scores.csv',uncertainty['source_score_hashes']['confirmation'])
    assert corrected['settings_unchanged'] and corrected['count_depths_exactly_match_actual_templates']
    packages=['numpy','scipy','pandas','h5py','anndata','scikit-learn','threadpoolctl','zstandard','vcc-cli','matplotlib']
    build.write_json(A/'runtime_versions.json',{'python':sys.version,'platform':platform.platform(),
        'packages':{name:version(name) for name in packages},
        'scope':'Runtime used for the final provenance audit and generated submission.'})
    build.write_json(A/'linked_provenance_validation.json',{'passed':True,'checked_files':sorted(set(checked)),
        'script_sha256':build.sha256(ROOT/'audit_v10_provenance.py'),'package_sha256':ready['sha256'],
        'runtime_versions_sha256':build.sha256(A/'runtime_versions.json'),
        'all_packaging_depth_corrected_diagnostics_and_frozen_code_links_match':True,
        'scope':'Record integrity after independent full-payload validation; no tuning or hidden-performance claim.'})
    print('V10 final linked provenance passed:',len(set(checked)),'files',flush=True)

if __name__=='__main__':main()
