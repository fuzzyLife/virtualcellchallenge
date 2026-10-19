"""Seal/check V69 root code, compact evidence, arrays and release artifacts."""
import argparse,json
from pathlib import Path
from v57_pie import read,write,sha
import v69_count_bridge as m


def paths():
    root=sorted(p for p in m.ROOT.glob('*v69*.py') if p.is_file())
    docs=[m.ROOT/'V69_STRATEGY.md',m.ROOT/'V69_RESULTS.md']
    evidence=sorted(p for p in m.A.rglob('*') if p.is_file() and p.name!='release_manifest.json')
    return root+docs+evidence+[m.OUT/'prediction.h5ad',m.OUT/'prediction.vcc']


def main(check=False):
    destination=m.A/'release_manifest.json'
    if check:
        manifest=read(destination)
        for i,r in enumerate(manifest['files']):
            p=m.ROOT/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],str(p)
            if (i+1)%50==0:print('V69 seal verified',i+1,'/',len(manifest['files']),flush=True)
        assert {p.relative_to(m.ROOT).as_posix() for p in paths()}=={r['path'] for r in manifest['files']}
        print('V69 release seal passed',len(manifest['files']),'files',flush=True);return
    assert not destination.exists();assert read(m.A/'independent_validation.json')['passed'] and read(m.A/'model_audit.json')['passed']
    records=[]
    for i,p in enumerate(paths()):
        records.append(dict(path=p.relative_to(m.ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
        if (i+1)%50==0:print('V69 seal hashed',i+1,flush=True)
    write(destination,dict(files=records,source_model_unchanged='V66 means; additional count bridge',uploaded=False,private_score=None,README_excluded_mutable_status=True))
    print('V69 release seal created',len(records),'files',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
