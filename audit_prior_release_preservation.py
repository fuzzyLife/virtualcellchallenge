"""Check immutable prior-release inventory and metadata without rehashing large arrays."""
from pathlib import Path
from v57_pie import read,write,sha
from seal_v69_release import paths
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/post_v69_review/analysis'
def main():
 A.mkdir(parents=True,exist_ok=True);manifest=read(ROOT/'experiments/v69_count_bridge/analysis/release_manifest.json');expected={r['path'] for r in manifest['files']};actual={p.relative_to(ROOT).as_posix() for p in paths()};assert expected==actual
 checked=[]
 for r in manifest['files']:
  p=ROOT/r['path'];assert p.stat().st_size==r['bytes']
  if p.suffix.lower() in ['.py','.md','.json','.csv','.txt']:
   assert sha(p)==r['sha256'],str(p);checked.append(r['path'])
 result=dict(passed=True,V69_inventory_files=len(expected),all_recorded_sizes_match=True,metadata_code_and_text_hashes_checked=len(checked),large_binary_hashes_rechecked=False,scope='Inventory, all sizes, and code/text/metadata hashes only. This is not a new full binary scientific-seal audit.',auditor_sha256=sha(__file__))
 write(A/'prior_release_preservation.json',result);print(result,flush=True)
if __name__=='__main__':main()
