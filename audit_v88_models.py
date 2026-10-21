"""Independent V88 convergence, exact likelihood and excluded-study source audits."""
from unittest.mock import patch
import argparse,gc
import numpy as np
from scipy import stats
from scipy.optimize import brentq
from scipy.special import ndtr
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import read,write,sha

def main(production=False):
 m.check();d=m.data();records=[]
 for outer in (['production'] if production else ['H1','Orion','CD4']):
  b=m.old.bank(outer)
  for fake in ([False] if production else [False,True]):
   net=m.model(outer,fake);r=read(m.A/'models'/(outer+('_shuffled' if fake else '')+'.json'));y,v,ok=m.old.observations(b,fake);ids=net['sample_indices'];ll,_,_=m.likelihood(y[ids],v[ids],ok[ids],b['similarity']);ix=np.array([i for i,(s,l) in enumerate(m.COMPONENTS) if s==0 or l==999]);details=[]
   for weights,axis in [(net['weights'],np.arange(len(m.COMPONENTS))),(net['shared_weights'][ix],ix)]:
    a=np.exp(ll[:,axis]-ll[:,axis].max(1,keepdims=True));g=a.T@(1/(a@weights))+.5/weights;lam=len(a)+.5*len(weights);res=float(abs(g-lam).max()/lam);assert res<1e-8;details.append(res)
   small=ids[:24];actual,_,_=m.likelihood(y[small],v[small],ok[small],b['similarity']);expected,_,_=m.old.likelihood(y[small],v[small],ok[small],b['similarity']);np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-8)
   for ci,(scale,length) in enumerate(m.COMPONENTS):
    for j,row in enumerate(small[:6]):
     use=np.flatnonzero(ok[row]);cov=np.diag(v[row,use].astype(float))
     if scale:cov+=scale**2*m.old.kernel(b['similarity'],length)[np.ix_(use,use)]
     val=stats.multivariate_normal.logpdf(y[row,use],mean=np.zeros(len(use)),cov=cov);np.testing.assert_allclose(actual[j,ci],val,rtol=1e-9,atol=1e-7)
   # An actual query feature vector is held control expression only, never treatment.
   qp=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/({'H1':'H1','Orion':'HCT116','CD4':'CD4_Rest_D1'}.get(outer,'H1')+'.npz')
   with np.load(qp,allow_pickle=False) as z:q=z['query'].copy()
   axis=net['axis'];known=np.isfinite(q[axis]);feature=np.zeros(len(axis));rank=stats.rankdata(q[axis][known]);rank-=rank.mean();feature[known]=rank/np.linalg.norm(rank);sim=net['rank_features']@feature;ll0,mu,var=m.likelihood(y[small],v[small],ok[small],b['similarity'],sim);prob=m.probabilities(ll0,net['weights'])
   for level in [.1,.9]:
    got=m.quantile(mu,var,prob,level)
    for j in range(len(small)):
     positive=var[j]>0
     def continuous(x):return np.sum(prob[j,positive]*ndtr((x-mu[j,positive])/np.sqrt(var[j,positive])))
     lo=continuous(0);atom=prob[j,~positive].sum()
     if lo<level<=lo+atom:want=0.
     elif level<=lo:want=brentq(lambda x:continuous(x)-level,-100,0)
     else:want=brentq(lambda x:continuous(x)+atom-level,0,100)
     np.testing.assert_allclose(got[j],want,rtol=0,atol=3e-8)
   c=r['comparison_to_V86'];assert c is None or c['converged_shifted_objective']>=c['V86_shifted_objective']-1e-7
   records.append(dict(outer=outer,shuffled=fake,KKT_relative_residuals=details,likelihood_coordinates=24,actual_Brent_quantiles=48,source_objective_noninferior_to_V86=True));del y,v,ok;gc.collect()
  if not production:
   excluded=np.array([m.old.GROUP[str(f)]==outer for f in d['row_families']]);poison=dict(d)
   for k in ['truth','errors','valid','knockdown']:poison[k]=d[k].copy();poison[k][excluded]=True if k=='valid' else 999
   poison['profiles']=d['profiles'].copy();poison['profiles'][np.unique(d['row_contexts'][excluded])]=999999;m.old.bank.cache_clear()
   with patch.object(m.old,'data',return_value=poison):
    other=m.old.bank(outer)
    for k in ['y','variance','valid','source_profiles','rank_features','axis','similarity','source_indices','source_names']:np.testing.assert_array_equal(b[k],other[k])
   m.old.bank.cache_clear();gc.collect()
  print('V88 source-model audit',outer,flush=True)
 write(m.A/('production_model_audit.json' if production else 'model_audit.json'),dict(passed=True,records=records,all_weights_converged=True,actual_mixture_quantiles_checked=True,whole_recipient_study_poison_has_no_effect=not production,script_sha256=sha(__file__)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--production',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):main(a.production)
