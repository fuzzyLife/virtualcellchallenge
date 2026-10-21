"""V78: calibrate cell variance while preserving exact V66 count margins."""
import argparse,gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v69_count_bridge as bridge
from v29_methods import project
from v57_pie import read,write,sha
from v62_conserved_response import data
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
ROOT=bridge.ROOT;OUT=ROOT/'experiments/v78_margin_variance';A=OUT/'analysis'
ARMS=('control10','control25','learned25','query50','query100')

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v78_margin_variance.py','v69_count_bridge.py','v29_methods.py','v29_pairwise_repair.py','v29_swap_kernel.py','v26_methods.py','v16_methods.py']
 inputs=['experiments/v66_assay_strength/analysis/independent_validation.json','experiments/v66_assay_strength/analysis/reference_model_arrays.npz']
 for name in bridge.SOURCES:inputs.extend([f'experiments/v69_count_bridge/analysis/source_shapes/{name}.npz',f'experiments/v69_count_bridge/analysis/source_shapes/{name}.json'])
 write(A/'plan.json',dict(version=78,dependencies={n:sha(ROOT/n) for n in deps},inputs={n:sha(ROOT/n) for n in inputs},arms=ARMS,hypothesis='V69 changed means and regressed. V70-V77 response mean corrections fail discrimination. Keep every realized V66 row total and every gene total EXACTLY, learn only cell variance. NB moment-matched shape bridge at zero mean effect, followed by reviewed V29 exact-margin projection and CPM repair. Related to old fixed-mean variation tests, but protects actual published V66 profiles and uses a count-shape proposal rather than latent neural residuals.',shape='Fit depth-conditioned NB alpha on V66 population and full query NTC controls. Reliable nonpanel genes require >=80 molecules in both. Restore query alpha; learned arm adds quarter of study-balanced V69 source dispersion shift with whole-study exclusion. Shape log-ratio clipped +/-ln4. Mix proposed counts with original at10%/25% before projection. Other genes receive zero shape change. Project to exact original row and column counts; retain original if global CPM drift>0.0005 or movement>15% of original molecules. Undefined alpha uses original; no target treatment input.',controls='Each arm uses the same shaping/rounding seed. Control arms compare alpha shuffled among reliable genes within ten expression bins. Learned arm compares cyclic target-shuffled V69 dispersion evidence at same strength. Zero shape is exact count identity. Existing V29 projector unchanged, including5% CPM repair budget.',development='First TWO unique fixed cases/context (10 populations, PDS>=2). Three arms versus V66 and respective shuffled controls. Profile/PDS must be bitwise identical, CPM bound passes. Select highest fidelity among candidates with fidelity gain>=.002, Jaccard gain>=.0005, reach no worse .005, NMAE available-study error no worse1%, and fidelity/Jaccard better than its shuffled control. Missing studies explicitly remain missing; available NMAE is a partial guard, never full six-metric support.',confirmation='Remaining unique cases after first2 (2/2/4/6/6,20 populations), two seeds0/+620000. No promotion after confirmation. Critical release requires exact raw profiles/PDS, mean CPM bound, fidelity/Jaccard favorable versus V66 and shuffle, reach no material regression, available-study NMAE<=1.01. Corrected expression numerator/jackknife may change and is NOT protected by bulk identity. All public biological data reused; no hidden-score forecast.',D3_opened=False,automatic_upload=False))
 print('V78 exact-margin variance protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/n)==h,n
 return p

def permuted_alpha(alpha,p,eligible):
 result=alpha.copy();ix=np.flatnonzero(eligible);order=ix[np.argsort(p[ix],kind='stable')]
 for block in np.array_split(order,10):
  if len(block)>1:result[block]=alpha[np.roll(block,1)]
 return result

def shape_proposal(parent,params,target_alpha,seed):
    """Independent mean-preserving Gamma coupling; no fourfold cap on query targets."""
    x=parent.toarray().astype(np.int64);result=x.copy();rng=np.random.default_rng(seed)
    depth=x.sum(1);p=params['p'];r0=1/params['alpha'];r1=1/target_alpha
    for lo in range(0,x.shape[1],512):
        hi=min(lo+512,x.shape[1]);ix=np.arange(lo,hi);ix=ix[(p[ix]>0)&(abs(r1[ix]-r0[ix])>1e-12)]
        if not len(ix):continue
        obs=x[:,ix];mu=depth[:,None]*p[ix];a=r0[ix][None,:];b=r1[ix][None,:]
        latent=rng.gamma(a+obs,mu/(a+mu));shifted=latent.copy();higher=r1[ix]>r0[ix];lower=r1[ix]<r0[ix]
        if higher.any():shifted[:,higher]+=rng.gamma((b-a)[:,higher],mu[:,higher]/a[:,higher])
        if lower.any():shifted[:,lower]*=rng.beta(b[:,lower],(a-b)[:,lower],size=(len(x),int(lower.sum())))
        shifted*=a/b;keep=np.minimum(np.divide(shifted,latent,out=np.ones_like(latent),where=latent>0),1)
        result[:,ix]=obs-rng.binomial(obs,np.clip(1-keep,0,1))+rng.poisson(np.maximum(shifted-latent,0))
    return sparse.csr_matrix(result)

def transform(parent,query,delta,arm,seed,shuffle=False,panel=None):
    parent=sparse.csr_matrix(parent,dtype=np.int32);original=parent.toarray();params=bridge.parameters(parent)
    eligible=(query['total']>=80)&(params['total']>=80)
    if panel is not None:eligible&=~panel
    target=np.asarray(query['alpha'],float)
    if shuffle and (arm.startswith('control') or arm.startswith('query')):target=permuted_alpha(target,np.asarray(query['p']),eligible)
    shift=np.log(target/params['alpha'])
    if arm=='learned25':shift=shift+.25*np.asarray(delta,float)
    limit=np.log(4000) if arm.startswith('query') else np.log(4)
    shift=np.where(eligible,np.clip(shift,-limit,limit),0)
    if not np.any(shift):return parent.copy(),dict(zero_shape_identity=True,bulk_exact=True,row_libraries_exact=True,CPM_L1=0.,fallback=False,changed_fraction=0.)
    if arm.startswith('query'):
        target_alpha=np.clip(params['alpha']*np.exp(shift),.005,20)
        proposed=shape_proposal(parent,params,target_alpha,seed)
    else:proposed,_=bridge.bridge(parent,np.zeros(parent.shape[1]),None,seed,params,shift,strength=1.)
    fraction={'control10':.1,'control25':.25,'learned25':.25,'query50':.5,'query100':1.}[arm]
    initial=(1-fraction)*original+fraction*proposed.toarray()
    count,diag=project(initial,parent,seed+780078)
    np.testing.assert_array_equal(np.asarray(count.sum(0)).ravel(),np.asarray(parent.sum(0)).ravel())
    np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),np.asarray(parent.sum(1)).ravel())
    moved=float(abs(count-parent).sum()/2/parent.sum());fallback=diag['mean_reference_L1']>.0005 or moved>(.35 if arm.startswith('query') else .15)
    if fallback:count=parent.copy()
    libs=np.asarray(parent.sum(1)).ravel();mu0=np.asarray(parent.multiply((1/libs)[:,None]).mean(0)).ravel();mu1=np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()
    achieved=bridge.parameters(count);mask=eligible
    before=float(abs(np.log(params['alpha'][mask]/query['alpha'][mask])).mean()) if mask.any() else None
    after=float(abs(np.log(achieved['alpha'][mask]/query['alpha'][mask])).mean()) if mask.any() else None
    diagnostics=dict(diag,alpha_log_error_before=before,alpha_log_error_after=after,zero_shape_identity=False,bulk_exact=True,row_libraries_exact=True,CPM_L1=float(abs(mu1-mu0).sum()),fallback=bool(fallback),changed_fraction=float(abs(count-parent).sum()/2/parent.sum()),reliable_genes=int(eligible.sum()))
    return count,diagnostics


def query(name):
 p=A/'query_params'/(name+'.npz');r=read(p.with_suffix('.json'));assert sha(p)==r['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def prepare_queries():
 check();from v25_methods import raw
 d=data();genes=d['genes'];pdir=A/'query_params';pdir.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family'];p=pdir/(name+'.npz')
  if p.exists():assert sha(p)==read(p.with_suffix('.json'))['sha256'];continue
  x,meta,measured,source=raw(name,genes);ids=np.asarray(ctx['left']);assert (meta.gene_target.to_numpy(str)[ids]=='Non-Targeting').all()
  fit=bridge.parameters(x[ids]);np.savez_compressed(p,alpha=fit['alpha'],p=fit['p'],total=fit['total'],genes=genes)
  write(p.with_suffix('.json'),dict(sha256=sha(p),source=name,controls=len(ids),source_record=source,script_sha256=sha(__file__)))
  print('V78 query variance',name,len(ids),flush=True);del x;gc.collect()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  else:prepare_queries()
