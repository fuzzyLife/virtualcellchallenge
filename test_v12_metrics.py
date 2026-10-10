"""Independent official DE metric and actual-count CPU Wilcoxon parity checks."""
import json,sys,subprocess,platform
from pathlib import Path
import numpy as np
from scipy import sparse
import build_submission as build
from v12_cpu_metrics import de,direction_metrics
from test_evaluation_criteria import official_modules,fixture,evaluate
A=build.ROOT/'experiments/v12_distribution_calibrated/analysis'
def main():
 A.mkdir(exist_ok=True,parents=True);dm,direction=official_modules();cases=[]
 for n,real,pred in [(9,9,9),(10,10,10),(14,4,14),(14,14,4),(14,0,0),(14,14,14)]:
  genes,truth,guess=fixture(n,real,pred)
  if (n,real,pred)==(14,14,14):truth['lfc'][:]=0
  exact=evaluate(dm,direction,genes,truth,guess);local=direction_metrics(guess,truth,genes,0)
  for k,v in exact.items():np.testing.assert_allclose(local[k],v,equal_nan=True)
  cases.append({'genes':n,'real_calls':real,'pred_calls':pred,'all_four_match':True})
 rng=np.random.default_rng(2026100911);genes=np.array(['TARGET']+[f'G{i:04}' for i in range(127)])
 # Sparse ties, unequal depths, low-control genes, constant and all-zero genes.
 rate=np.exp(rng.uniform(-5,4,len(genes)));rate[1]=0;rate[2]=.00001
 c=rng.poisson(rng.lognormal(0,.3,(183,1))*rate);t=rng.poisson(rng.lognormal(0,.3,(117,1))*rate*np.r_[np.repeat(2.,30),np.ones(98)])
 c[:,3]=t[:,3]=10
 np.savez(A/'cpu_parity_fixture.npz',genes=genes,control=c,treated=t)
 local=de(sparse.csr_matrix(t),sparse.csr_matrix(c),genes)
 np.savez(A/'cpu_parity_local.npz',**local)
 cmd=[str(build.ROOT/'.v12_eval_env/Scripts/python.exe'),str(build.ROOT/'verify_v12_official_cpu.py')]
 subprocess.run(cmd,check=True)
 official=json.loads((A/'official_cpu_parity.json').read_text())
 assert official['passed']
 build.write_json(A/'metric_tests.json',{'passed':True,'fixtures':cases,'official_counts_pipeline':official,
  'local_python':platform.python_version(),'script_sha256':build.sha256(Path(__file__)),
  'methods_sha256':build.sha256(build.ROOT/'v12_cpu_metrics.py')})
 print('Six metric edge fixtures and official CPU actual-count pipeline parity passed',flush=True)
if __name__=='__main__':main()
