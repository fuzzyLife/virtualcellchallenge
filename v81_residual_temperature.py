"""V81: exact-margin count-residual contraction and opposite-direction controls."""
import argparse
import numpy as np
from scipy import sparse
import v78_margin_variance as old
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v81_residual_temperature';A=OUT/'analysis';ARMS=('contract25','contract50');NEG_SUFFIX='_opposite'
data=old.data;contexts=old.contexts;key=old.key;checked_folder=old.checked_folder;bridge=old.bridge;query=old.query

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v81_residual_temperature.py','v78_margin_variance.py','v29_methods.py','v29_pairwise_repair.py','v26_methods.py','v16_methods.py']
 write(A/'plan.json',dict(version=81,dependencies={n:sha(ROOT/n) for n in deps},arms=ARMS,hypothesis='V78 and partial V79/V80 results show restoring control variance can remove correct DE calls. Test the opposite decoder direction: contract V66 residuals toward the library-scaled bulk expectation, then enforce exact parent row/gene sums and CPM mean. Related to earlier dispersion tuning, but controls actual published V66 margins. This is a noise-calibration diagnostic, not a claim that low biological variance is realistic.',proposal='Initial=(1-gamma)*X+gamma*library*bulk; gamma .25/.50. Negative control uses gamma -.25/-.50 and clips negative proposals before the same projection. No treatment values fitted; no learned source coefficients. All genes receive the same global temperature. Fixed count sums protect raw bulk/PDS; uncertainty-corrected MSE may change.',selection='First2 unique cases/context. Exact profiles/PDS and CPM<=.0005; FID gain>=.002, JAC>=baseline-.0005, REACH>=baseline-.005, available-study NMAE<=1.01; genuine FID/JAC better than opposite noise expansion. Highest FID, tie JAC. This separately specified opposite-direction hypothesis does not promote failed V78/V79/V80 candidates.',confirmation='Remaining20 cases both seeds0/620000; no reselection, same guards. A zero-response null test must additionally show no material increase in false DE calls before any production release. All reused public biology; no hidden overall forecast.',automatic_upload=False))

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 old.check();return p

def transform(parent,query,delta,arm,seed,shuffle=False,panel=None):
 parent=sparse.csr_matrix(parent,dtype=np.int32);x=parent.toarray();libs=x.sum(1,dtype=np.int64);bulk=x.sum(0)/libs.sum();g={'contract25':.25,'contract50':.5}[arm];g=-g if shuffle else g
 initial=np.maximum((1-g)*x+g*libs[:,None]*bulk[None,:],0)
 count,diag=old.project(initial,parent,seed+810081)
 moved=float(abs(count-parent).sum()/2/parent.sum());fallback=diag['mean_reference_L1']>.0005 or moved>.35
 if fallback:count=parent.copy()
 np.testing.assert_array_equal(count.sum(0),parent.sum(0));np.testing.assert_array_equal(count.sum(1),parent.sum(1));p0=np.asarray(parent.multiply((1/libs)[:,None]).mean(0)).ravel();p1=np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()
 return count,dict(diag,CPM_L1=float(abs(p1-p0).sum()),fallback=bool(fallback),changed_fraction=float(abs(count-parent).sum()/2/parent.sum()),bulk_exact=True,row_libraries_exact=True,temperature_gamma=g)

if __name__=='__main__':freeze()
