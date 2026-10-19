"""Finalize assessment metadata after fixed V69 counts and confirmation complete."""
from pathlib import Path
import shutil
from v57_pie import read,write,sha
import v69_count_bridge as m


def main():
    path=m.A/'reference_model.json';record=read(path);assessment=read(m.A/'confirmation_assessment.json');matched=read(m.A/'matched_control_assessment.json');selection=read(m.A/'selection.json')
    assert (m.OUT/'prediction.h5ad').exists() and (m.A/'run_parameters.json').exists()
    assert assessment['selected']==record['candidate']==selection['selected']==matched['candidate']
    if not record.get('confirmation_pending',False):
        print('V69 confirmation metadata already finalized',flush=True);return
    arrays_sha=sha(m.A/'reference_model_arrays.npz');assert arrays_sha==record['arrays_sha256']
    runpath=m.A/'run_parameters.json';run=read(runpath);old_hash=sha(path);assert run['model_sha256']==old_hash
    folder=m.A/'execution_history';folder.mkdir(exist_ok=True)
    for p in (path,runpath):
        dest=folder/('before_confirmation_status_'+p.name);assert not dest.exists();shutil.copyfile(p,dest)
    supported=bool(assessment['release_supported'] and matched['beats_matched_control_profile'] and record['full_bridge_null_passed'] and record['candidate']!='hybrid25' and not selection['experimental_fallback'])
    record.update(confirmation_pending=False,status='public_supported' if supported else 'experimental_failed_public_guards',available_public_guards_passed=supported,assessment_sha256=sha(m.A/'confirmation_assessment.json'),matched_control_assessment_sha256=sha(m.A/'matched_control_assessment.json'),pre_confirmation_metadata_sha256=old_hash,scientific_arrays_unchanged_at_finalization=True)
    write(path,record);run.update(model_sha256=sha(path),pre_confirmation_model_metadata_sha256=old_hash,assessment_status_only_updated=True);write(runpath,run)
    assert sha(m.A/'reference_model_arrays.npz')==arrays_sha
    write(m.A/'confirmation_finalization.json',dict(scientific_arrays_sha256=arrays_sha,scientific_array_bytes_unchanged=True,counts_already_generated=True,previous_model_metadata_sha256=old_hash,final_model_metadata_sha256=sha(path),finalizer_sha256=sha(Path(__file__)),candidate_unchanged=True))
    print('V69 confirmation status finalized; model/count arrays unchanged',flush=True)

if __name__=='__main__':main()
