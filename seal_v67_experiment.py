"""Seal or verify the completed V67 experiment without altering prior releases."""
from pathlib import Path
import argparse
import v67_intervention_network as m
from v57_pie import read,write,sha


def main(check=False):
    m.check();path=m.A/'experiment_manifest.json'
    if check:
        receipt=read(path)
        for name,h in receipt['files_sha256'].items():assert sha(m.ROOT/name)==h,name
        audit=read(m.A/'results_audit.json')
        for name,h in audit['source_sha256'].items():assert sha(m.ROOT/name)==h,name
        for p in (m.A/'predictions').glob('*.json'):
            assert sha(p.with_suffix('.npz'))==read(p)['sha256'],str(p)
        assert read(m.A/'mean_results.json')['release_supported']==receipt['release_supported']
        print('V67 experiment seal verified',len(receipt['files_sha256']),'compact files and all prediction/audit input hashes.',flush=True)
        return
    assert not path.exists()
    result=read(m.A/'mean_results.json');assert read(m.A/'results_audit.json')['passed']
    chosen=set()
    for pattern in ('*v67*.py','V67*.md'):chosen.update(m.ROOT.glob(pattern))
    for p in m.A.rglob('*'):
        if p.is_file() and p.suffix in ('.json','.csv','.md','.png') and p!=path:chosen.add(p)
    write(path,dict(files_sha256={str(p.relative_to(m.ROOT)):sha(p) for p in sorted(chosen)},
          release_supported=result['release_supported'],private_score=None,
          submission_generated=False,uploaded=False,D3_opened=False,
          scope='Frozen V67 method, protocol, public response predictions and scoring receipts, negative-control comparisons and documentation. Prediction arrays are bound through sidecar and audit hashes; large arrays are ignored. README is mutable and excluded. Original V63 and V65/V66 seals are unchanged and checked separately.'))
    print('Sealed V67',len(chosen),'compact files.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
