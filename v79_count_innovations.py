"""V79: underdispersion-aware integer innovations, exact V66 profile margins."""
import argparse,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v78_margin_variance as old
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v79_count_innovations';A=OUT/'analysis'
ARMS=('nb_query','turnover50','turnover80','source50')
data=old.data;contexts=old.contexts;key=old.key;checked_folder=old.checked_folder;bridge=old.bridge;query=old.query

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v79_count_innovations.py','v78_margin_variance.py','v29_methods.py','v29_pairwise_repair.py','v29_swap_kernel.py','v26_methods.py','v16_methods.py']
 write(A/'plan.json',dict(version=79,dependencies={n:sha(ROOT/n) for n in deps},arms=ARMS,hypothesis='Code inspection: V66 initial count mixture is25%expectation+65%Poisson+10%empirical, so biological variance and even Poisson variance can be attenuated. NB shape floors cannot represent underdispersion. Use integer thinning plus NB innovations matching a target control-based variance, then preserve exact V66 count margins. Specified before inspecting any V78 treatment DE endpoint.',formula='For mean mu, parent Var=mu+a0*mu^2 with signed a0 (no positive floor), keep rho<=.8. Y=Binomial(X,rho)+NB((1-rho)*mu,a_add), a_add=(a_target-rho^2*a0)/(1-rho)^2. Adapt rho downward if needed so a_add>=0. Then VarY=mu+a_target*mu^2 under working depth-conditioned moments, for ANY parent dispersion. This is a moment model, not an inferred biological turnover rate or guaranteed marginal distribution after projection.',arms_description='NB query:rho0, full observed NTC alpha. turnover50/80:retain.5/.8 of parent molecules plus innovations. source50 adds quarter study-excluded V69 shape delta. NTC>=80/parent>=20 molecules; exclude only current own gene from direct proposal (other panel genes may be calibrated). Same query alpha shuffle within expression bins, or target-shuffled source evidence for source50.',constraints='Original V29 exact-row/exact-column/CPM projection unchanged. Revert to original parent if global CPM L1>.0005 or moved fraction>.35. Exact raw bulk/PDS invariant; corrected expression sampling terms and DE means are not invariant and must be reported. Projection can alter carrier genes.',selection='First2 unique fixed cases/context (10). All profile/PDS invariants must pass. Precision-focused criterion: fidelity gain>=.01, Jaccard gain>=.001, reach decrease<=.02, available-study NMAE<=1.01; genuine fidelity>shuffled and Jaccard>=shuffled. This is a distinct prospectively specified precision/overlap objective, not promotion of a failed V78 candidate. Choose maximum fidelity, tie Jaccard. Missing studies remain missing; no full-six-metric support claim.',confirmation='Remaining20 cases, both seeds0/+620000; no strength/model reselection. Require same precision/overlap/recall/NMAE and shuffle guards, exact row/gene/profile invariants and independent DE/rank checks before production. Public biology reused; no hidden overall-score estimate.',D3_opened=False,automatic_upload=False))
 print('V79 underdispersion-aware innovation protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 old.check();return p

def signed_alpha(parent):
 x=sparse.csr_matrix(parent,dtype=float);n=x.shape[0];depth=np.asarray(x.sum(1)).ravel();tot=np.asarray(x.sum(0)).ravel();p=tot/depth.sum();den=p*p*(depth@depth)
 residual=np.maximum(np.asarray(x.power(2).sum(0)).ravel()-2*p*np.asarray(x.T@depth).ravel()+den,0)*n/(n-1)
 return np.divide(residual-tot,den,out=np.zeros_like(p),where=den>0)

def innovation_parameters(a0,target,keep):
 rho=np.full_like(a0,keep,dtype=float);positive=a0>0
 rho[positive]=np.minimum(rho[positive],np.sqrt(np.maximum(target[positive],0)/a0[positive]))
 add=np.maximum((target-rho*rho*a0)/(1-rho)**2,0)
 return rho,add

def draw_innovations(obs,mu,rho,add,rng):
 kept=rng.binomial(obs,rho[None,:]);innovation_mean=(1-rho)[None,:]*mu;rate=innovation_mean.copy();positive=add>1e-10
 if positive.any():rate[:,positive]=rng.gamma(1/add[positive],innovation_mean[:,positive]*add[positive])
 return kept+rng.poisson(rate)

def transform(parent,query,delta,arm,seed,shuffle=False,panel=None):
 parent=sparse.csr_matrix(parent,dtype=np.int32);original=parent.toarray();fit=bridge.parameters(parent);a0=signed_alpha(parent)
 eligible=(query['total']>=80)&(fit['total']>=20)
 if panel is not None:eligible&=~panel
 if not eligible.any():return parent.copy(),dict(CPM_L1=0.,fallback=False,changed_fraction=0.,bulk_exact=True,row_libraries_exact=True,zero_shape_identity=True)
 target=np.asarray(query['alpha'],float).copy()
 if shuffle and arm!='source50':target=old.permuted_alpha(target,np.asarray(query['p']),eligible)
 if arm=='source50':target=np.clip(target*np.exp(.25*np.asarray(delta,float)),.005,20)
 keep={'nb_query':0.,'turnover50':.5,'turnover80':.8,'source50':.5}[arm]
 rho,add=innovation_parameters(a0,target,keep);proposal=original.copy();rng=np.random.default_rng(seed);libs=original.sum(1,dtype=np.int64)
 for lo in range(0,original.shape[1],512):
  ix=np.arange(lo,min(lo+512,original.shape[1]));ix=ix[eligible[ix]]
  if not len(ix):continue
  mu=libs[:,None]*fit['p'][ix];proposal[:,ix]=draw_innovations(original[:,ix],mu,rho[ix],add[ix],rng)
 count,diag=old.project(proposal,parent,seed+790079)
 np.testing.assert_array_equal(count.sum(0),parent.sum(0));np.testing.assert_array_equal(count.sum(1),parent.sum(1))
 moved=float(abs(count-parent).sum()/2/parent.sum());fallback=diag['mean_reference_L1']>.0005 or moved>.35
 if fallback:count=parent.copy()
 before=float(abs(np.log(fit['alpha'][eligible]/query['alpha'][eligible])).mean());after_alpha=bridge.parameters(count)['alpha'];after=float(abs(np.log(after_alpha[eligible]/query['alpha'][eligible])).mean())
 cpm0=np.asarray(parent.multiply((1/libs)[:,None]).mean(0)).ravel();cpm1=np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()
 return count,dict(diag,CPM_L1=float(abs(cpm1-cpm0).sum()),fallback=bool(fallback),changed_fraction=float(abs(count-parent).sum()/2/parent.sum()),bulk_exact=True,row_libraries_exact=True,zero_shape_identity=False,alpha_log_error_before=before,alpha_log_error_after=after,reliable_genes=int(eligible.sum()),signed_underdispersed_gene_fraction=float(np.mean(a0[eligible]<0)),nominal_retention=keep)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze']);a=p.parse_args();freeze()
