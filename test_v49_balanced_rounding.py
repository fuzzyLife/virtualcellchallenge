"""Compare bounded rounding with the previous greedy integer margin repair."""
import numpy as np
import v49_balanced_rounding as m
from v26_methods import integerize as greedy
from pathlib import Path
import build_submission as b
def main():
 rng=np.random.default_rng(49049);records=[];found=False
 for trial in range(30):
  n=20;g=30;libs=rng.integers(80,200,n);ref=np.array([rng.multinomial(int(l),np.ones(g)/g) for l in libs]);init=rng.gamma(2,2,(n,g));totals=ref.sum(0);e,steps=m.fit_expected(init,libs,totals)
  out,moves=m.integerize(e,libs,totals,trial);old,_=greedy(e/libs[:,None],libs,totals,trial);error=float(abs(out-e).max());old_error=float(abs(old-e).max());found|=old_error>1.01
  records.append({'trial':trial,'new_max_entry_error':error,'old_max_entry_error':old_error,'binary_transfers':moves,'IPF_steps':steps})
 assert found,'No greedy outlier encountered in comparison'
 A=b.ROOT/'experiments/v49_balanced_population/analysis';A.mkdir(parents=True,exist_ok=True);b.write_json(A/'rounding_contract_tests.json',{'passed':True,'trials':records,'all_rows_columns_exact':True,'every_new_entry_floor_or_ceil':True,'old_greedy_exceeds_one_count_error':found,'script_sha256':b.sha256(Path(__file__)),'rounder_sha256':b.sha256(Path(m.__file__))})
 print('V49 bounded-rounding tests PASS; greedy max error',max(r['old_max_entry_error'] for r in records),'new',max(r['new_max_entry_error'] for r in records),flush=True)
if __name__=='__main__':main()


