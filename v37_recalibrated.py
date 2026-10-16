"""Post-mixture dual-moment calibration; fixed corrective hypothesis, no weight grid."""
import numpy as np
from scipy import sparse
from v37_ensemble import emit as original
NAME='recalibrated_population50'

def emit(template,effect,base,retrieval,own,seed,prepared,atlas,arm):
 if arm=='v24':return original(template,effect,base,retrieval,own,seed,prepared,atlas,arm)
 assert arm==NAME
 mixed,prior=original(template,effect,base,retrieval,own,seed,prepared,atlas,'population50')
 libs,p,bulk=prepared;libs=np.asarray(libs,np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum()
 initial=mixed.toarray().astype(float)/libs[:,None]
 count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
 np.testing.assert_array_equal(count.sum(1),libs)
 rp=(count/libs[:,None]).mean(0);rb=count.sum(0)/libs.sum()
 return sparse.csr_matrix(count),{'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'precalibration_mean_L1':prior['mean_probability_L1'],'precalibration_bulk_L1':prior['bulk_probability_L1']}
