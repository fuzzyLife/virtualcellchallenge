"""An unseen value in a train-constant donor slot must not activate random weights."""
import numpy as np
import v52_joint_response as old
from run_v53_joint_guard import normalizer,transform

def main():
 train=np.array([[1.,np.nan,7.],[2.,np.nan,7.],[3.,np.nan,7.]],np.float32)
 query=np.array([[2.,.4,99.]],np.float32);med,scale=normalizer(train);out=transform(query,med,scale)
 assert np.isinf(scale[1:]).all();np.testing.assert_array_equal(out[0,[1,2,4,5]],0)
 changed=query.copy();changed[0,1]=-.5;changed[0,2]=-99;np.testing.assert_array_equal(transform(changed,med,scale),out)
 before_med,before_scale=old.normalizer(train);before=old.transform(query,before_med,before_scale)
 assert before[0,1]!=0 and before[0,2]!=0
 np.testing.assert_allclose(out[:,0],before[:,0]);assert np.isfinite(out).all()
 # Active-feature missingness remains explicit.
 train=np.array([[1.,1.],[2.,np.nan],[3.,4.]],np.float32);med,scale=normalizer(train);q=transform(np.array([[np.nan,2.]]),med,scale);assert q[0,2]==1
 print('V53 constant/unseen donor-channel, active-feature parity and missingness checks passed',flush=True)
if __name__=='__main__':main()
