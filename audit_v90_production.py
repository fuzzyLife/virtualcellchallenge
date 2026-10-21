"""Check V90 production finite-iteration prior, without claiming optimizer convergence."""
import numpy as np
from scipy import stats
from threadpoolctl import threadpool_limits
import v90_diagnostic_release as m
from v57_pie import read,write,sha
if __name__=='__main__':
 with threadpool_limits(limits=2):
  m.check();impl=m.implementation();b=impl.bank('production');records=[]
  for fake in [False,True]:
   net=impl.model('production',fake);r=read(m.A/'models'/('production'+('_shuffled' if fake else '')+'.json'));y,v,ok=impl.observations(b,fake);ii=net['sample_indices'][:24];a,_,_=impl.likelihood(y[ii],v[ii],ok[ii],b['similarity']);gold,_,_=m.old.likelihood(y[ii],v[ii],ok[ii],b['similarity']);np.testing.assert_allclose(a,gold,rtol=1e-10,atol=1e-8)
   for weights in [net['weights'],net['shared_weights']]:assert (weights>=0).all();np.testing.assert_allclose(weights.sum(),1,atol=1e-12)
   for ci,(scale,length) in enumerate(impl.COMPONENTS):
    for j,row in enumerate(ii[:5]):
     use=np.flatnonzero(ok[row]);cov=np.diag(v[row,use].astype(float))
     if scale:cov+=scale**2*impl.kernel(b['similarity'],length)[np.ix_(use,use)]
     val=stats.multivariate_normal.logpdf(y[row,use],cov=cov);np.testing.assert_allclose(a[j,ci],val,rtol=1e-9,atol=1e-7)
   records.append(dict(shuffled=fake,iterations=len(r['em_history']),final_weight_change=r['em_history'][-1]['maximum_weight_change'],generic_likelihood_parity=True,convergence_claim=False))
  write(m.A/'production_model_audit.json',dict(passed=True,records=records,source_families=b['source_names'].tolist(),all_five_measured_source_families_used=True,script_sha256=sha(__file__)));print('V90 production prior parity audit passed',flush=True)
