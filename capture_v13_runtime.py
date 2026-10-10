"""Record actual production dependencies without changing the installed environment."""
from pathlib import Path
import importlib.metadata as metadata,json,platform,sys
from v13_methods import ROOT,A
import build_submission as build

def main():
 names=['numpy','scipy','pandas','anndata','h5py','scikit-learn','threadpoolctl','requests','zstandard']
 mapping=metadata.packages_distributions();names+=mapping.get('vcc',[])
 versions={name:metadata.version(name) for name in sorted(set(names))}
 file=ROOT/'v13_requirements.txt';file.write_text(''.join(f'{k}=={v}\n' for k,v in versions.items()),encoding='utf-8')
 build.write_json(A/'runtime.json',{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':versions,
  'requirements_sha256':build.sha256(file),'official_comparator':'Isolated v12_official_requirements.txt environment and pinned cell-eval2 commit; public anchor CSVs reused only for identical actual panels.',
  'script_sha256':build.sha256(Path(__file__))})
 print(versions,flush=True)
if __name__=='__main__':main()
