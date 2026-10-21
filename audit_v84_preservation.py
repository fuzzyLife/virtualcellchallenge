"""Check the existing V69 release inventory without rehashing its large payloads."""
from pathlib import Path
from v57_pie import read,write,sha
from seal_v69_release import paths
import v84_intervention_operator as m
manifest=read(m.ROOT/'experiments/v69_count_bridge/analysis/release_manifest.json')
expected={r['path'] for r in manifest['files']};actual={p.relative_to(m.ROOT).as_posix() for p in paths()};assert expected==actual
checked=[]
for r in manifest['files']:
 p=m.ROOT/r['path'];assert p.stat().st_size==r['bytes']
 if p.suffix.lower() in ['.py','.md','.json','.csv','.txt']:
  assert sha(p)==r['sha256'],str(p);checked.append(r['path'])
write(m.A/'prior_release_preservation.json',dict(passed=True,V69_inventory_files=len(expected),all_recorded_sizes_match=True,metadata_code_text_hashes_checked=len(checked),large_binary_hashes_rechecked=False,scope='Inventory and all sizes, plus code/text/metadata hashes. Not a new full binary scientific-seal audit.',auditor_sha256=sha(__file__)))
print('Prior release inventory preserved',len(expected),'files;',len(checked),'metadata/code hashes checked',flush=True)
