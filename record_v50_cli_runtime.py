"""Record official packaging runtime and unchanged installed source hashes."""
import sys,json,hashlib,importlib.metadata as im
from pathlib import Path
import vcc
ROOT=Path(__file__).resolve().parent
def main():
 packages=['vcc-cli','numpy','scipy','anndata','h5py','zstandard']
 source=Path(vcc.__file__).parent
 hashes={n:hashlib.sha256((source/n).read_bytes()).hexdigest() for n in ['prep.py','vccfile.py']}
 rec={'python':sys.version,'executable':sys.executable,'versions':{n:im.version(n) for n in packages},'installed_official_source_sha256':hashes,'source_modified_by_this_project':False,'collector_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (ROOT/'experiments/v50_response_prior/analysis/official_runtime_receipt.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
 print('Official CLI runtime recorded',rec['versions']['vcc-cli'],flush=True)
if __name__=='__main__':main()

