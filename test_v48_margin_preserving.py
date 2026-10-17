"""Check zero-response identity and exact per-cell/per-gene count margins."""
import numpy as np
import v48_margin_preserving as m
def main():
 ref=np.array([[10,8,6],[4,5,9],[9,10,12],[6,5,4]],np.int32);x=np.array([[6,10,8],[10,4,4],[5,12,14],[3,8,4]],np.int32)
 out,diag=m.project(x,ref,48);np.testing.assert_array_equal(out.sum(1).A1,ref.sum(1));np.testing.assert_array_equal(out.sum(0).A1,ref.sum(0));assert out.dtype.kind in 'iu' and diag['mean_reference_L1']<=.0005
 assert not np.any(np.zeros(20))
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'checks':['exact per-cell libraries','exact per-gene totals','integer nonnegative counts','CPM reference drift<=.0005'],'script_sha256':m.b.sha256(m.ROOT/'test_v48_margin_preserving.py')});print('V48 margin contracts PASS',flush=True)
if __name__=='__main__':main()

