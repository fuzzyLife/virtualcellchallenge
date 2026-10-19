"""Record the exact existing scientific runtime used for V65/V66 computations."""
from pathlib import Path
import sys,platform,importlib.metadata
import numpy as np
from scipy import stats
from threadpoolctl import threadpool_info,threadpool_limits
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent

def main():
 d=read(ROOT/'experiments/v65_closed_response/analysis/release_selection.json');version=d['selected_version'];assert version in (65,66);folder='v65_closed_response' if version==65 else 'v66_assay_strength';A=ROOT/'experiments'/folder/'analysis'
 names=('numpy','scipy','pandas','anndata','h5py','threadpoolctl','torch','scikit-learn','matplotlib','zstandard');versions={}
 for name in names:
  try:versions[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:versions[name]=None
 with threadpool_limits(limits=2):libraries=threadpool_info()
 write(A/'runtime_versions.json',dict(python=sys.version,executable=sys.executable,prefix=sys.prefix,platform=platform.platform(),packages=versions,native_libraries_under_scientific_thread_limit=libraries,private_generation_workers=4,public_generation_seed_offsets=[0,620000],production_seed=20261011,official_cli_version_scope='Read exact CLI version/schema from submission_ready.json package_meta after official local prep.',recorder_sha256=sha(__file__)))
 print('Recorded selected V'+str(version)+' scientific runtime:',versions,flush=True)
if __name__=='__main__':main()
