"""Seal/verify V68's completed scientific test while preserving prior releases."""
import argparse
from pathlib import Path
import v68_population_balancing as m
from v57_pie import read,write,sha


def verify_arrays():
    for p in (m.A/'frames').glob('*/receipt.json'):
        r=read(p);assert sha(p.parent/'controls.npz')==r['controls_sha256'];assert sha(p.parent/'cells.csv')==r['cells_sha256']
    for dirname in ('counts','fits','null','primary_null'):
        for p in (m.A/dirname).glob('*.json'):
            r=read(p)
            if 'arrays_sha256' not in r:continue
            if r.get('arm')=='v63' and 'reused_source' in r:
                prior=m.ROOT/r['reused_source'];assert sha(prior)==r['reused_receipt_sha256'];array=prior.with_suffix('.npz')
            else:array=p.with_suffix('.npz')
            assert sha(array)==r['arrays_sha256'],str(p)
    extra=read(m.A/'primary_null_plan.json');assert extra['original_plan_sha256']==sha(m.A/'plan.json') and extra['helper_sha256']==sha(m.ROOT/'run_v68_primary_null.py')


def main(check=False,refresh_docs=False):
    m.check();verify_arrays();path=m.A/'experiment_manifest.json'
    audit=read(m.A/'results_audit.json');assert audit['passed']
    for name,h in audit['source_sha256'].items():assert sha(m.ROOT/name)==h,name
    if check:
        receipt=read(path)
        for name,h in receipt['files_sha256'].items():assert sha(m.ROOT/name)==h,name
        print('V68 scientific seal verified',len(receipt['files_sha256']),'compact files, frames, populations and audit input hashes.',flush=True);return
    previous_sha=None
    if refresh_docs:
        previous=read(path);previous_sha=sha(path)
        allowed={'V68_RESULTS.md','document_v68_results.py','seal_v68_experiment.py'}
        for name,h in previous['files_sha256'].items():
            if name not in allowed:assert sha(m.ROOT/name)==h,name
        archive=m.A/'documentation_history'/('revision_'+previous_sha[:12]+'_manifest.json')
        assert not archive.exists();archive.write_bytes(path.read_bytes())
        prior_revision=m.A/'documentation_revision.json'
        if prior_revision.exists():
            backup=archive.with_name('revision_'+previous_sha[:12]+'_metadata.json');assert not backup.exists();backup.write_bytes(prior_revision.read_bytes())
        write(m.A/'documentation_revision.json',dict(previous_manifest_sha256=previous_sha,allowed_documentation_files=sorted(allowed),scientific_inputs_and_outcomes_unchanged=True,reason='Clarify that the displayed profile errors are squared-error values; fix spacing. Preserve the initial report, renderer, sealer and manifest.'))
    else:assert not path.exists()
    decision=read(m.A/'final_decision.json');assert not decision['submission_generated']
    chosen=set()
    for pattern in ('*v68*.py','V68*.md'):chosen.update(m.ROOT.glob(pattern))
    for p in m.A.rglob('*'):
        if p.is_file() and p.suffix in ('.json','.csv','.md','.png','.py') and p!=path:chosen.add(p)
    write(path,dict(files_sha256={str(p.relative_to(m.ROOT)):sha(p) for p in sorted(chosen)},
          superseded_documentation_manifest_sha256=previous_sha,release_supported=decision['release_supported'],submission_generated=False,uploaded=False,D3_opened=False,
          scope='Frozen V68 population-sampling methods, protocol and stricter pre-null addendum, paired development and null receipts, independent count/rank audit, plots and reports. Large arrays bound through receipt hashes and ignored. README mutable and excluded. Prior releases verified separately.'))
    print('Sealed V68',len(chosen),'compact files.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--refresh-docs',action='store_true');a=p.parse_args();assert not(a.check and a.refresh_docs);main(a.check,a.refresh_docs)


