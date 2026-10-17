"""Response-identity calibrator: zero identity, bounds, sign and masked-label tests."""
import numpy as np,torch
from threadpoolctl import threadpool_limits
import v40_identity_calibration as m
def main():
 rng=np.random.default_rng(40);base=rng.normal(0,.2,(8,128)).astype(np.float32);base[:,0]=0;theta=rng.normal(0,3,(2,128)).astype(np.float32)
 np.testing.assert_array_equal(m.predict(theta,base,0),base)
 for alpha in [.1,.25]:
  out=m.predict(theta,base,alpha);assert (out*base>=0).all();np.testing.assert_array_equal(out[:,0],base[:,0]);ratio=out[:,1:]/base[:,1:];assert ratio.min()>=1-alpha-1e-6 and ratio.max()<=1+alpha+1e-6
 p=torch.nn.Parameter(torch.zeros((2,128)));mask=torch.ones((8,128),dtype=torch.bool);mask[:,:3]=False
 group={'base':torch.tensor(base),'truth':torch.tensor(base.copy()),'mask':mask,'axis':torch.ones(128,dtype=torch.bool),'query':torch.ones(128)*20}
 ll,ce,_=m.loss(p,group,np.arange(8));wrong=m.loss(p,group,np.arange(8),True)[1];assert float(ce)<float(wrong)
 ll.backward();assert torch.isfinite(p.grad).all();value=float(ll);grad=p.grad.detach().clone();p.grad=None;group['truth'][:,:3]=2.5
 other=m.loss(p,group,np.arange(8))[0];other.backward();assert float(other)==value;torch.testing.assert_close(p.grad,grad,rtol=0,atol=0)
 m.write(m.A/'contract_tests.json',{'passed':True,'zero_gain_bitwise_baseline_identity':True,'positive_bounded_sign_preserving_gain':True,'masked_labels_do_not_affect_loss_or_gradient':True,'correct_target_contrast_beats_wrong_target':True,'finite_autograd':True,'script_sha256':m.b.sha256(m.ROOT/'test_v40_identity_calibration.py')})
 print('V40 calibration contracts PASS',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()

