"""Record the actual scientific runtime for the V50 experiment and generator."""
import sys,platform,importlib.metadata as im
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_info
import build_submission as b
def main():
 packages=['numpy','scipy','pandas','anndata','h5py','scikit-learn','joblib','threadpoolctl','torch','numba','llvmlite']
 b.write_json(b.ROOT/'experiments/v50_response_prior/analysis/runtime_receipt.json',{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'versions':{p:im.version(p) for p in packages},'threadpools':threadpool_info(),'BLAS_threads_during_experiment':2,'requirements_sha256':b.sha256(b.ROOT/'requirements-v35-science.txt'),'script_sha256':b.sha256(Path(__file__))})
 print('V50 scientific runtime recorded',flush=True)
if __name__=='__main__':main()

