"""Restore exported V61/V59 caches after temporary-directory cleanup; no overwrite."""
import shutil,tempfile
from pathlib import Path
from v57_pie import ROOT,read,sha

def main():
 folder=ROOT/'experiments/v61_pooled_state_response/analysis';manifest=read(folder/'array_exports.json');restored=0
 for row in manifest['files']:
  rel=Path(row['file']);source=folder/'arrays'/rel;assert sha(source)==row['sha256']
  if rel.parts[0]=='reused_raw':base=Path(tempfile.gettempdir())/'vcc_v59_state_composition';dest=base/'raw'/rel.name
  else:base=Path(tempfile.gettempdir())/'vcc_v61_pooled_state_response';dest=base/rel
  assert dest.resolve().is_relative_to(base.resolve())
  if dest.exists():assert sha(dest)==row['sha256']
  else:dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest);assert sha(dest)==row['sha256'];restored+=1
 print('Verified',len(manifest['files']),'exported arrays; restored',restored,flush=True)
if __name__=='__main__':main()
