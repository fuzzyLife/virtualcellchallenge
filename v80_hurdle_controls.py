"""V80: empirical-zero calibrated hurdle counts, exact published V66 margins."""
import argparse,gc,json
import numpy as np
from scipy import sparse
from scipy.special import expit,logit
from threadpoolctl import threadpool_limits
import v78_margin_variance as old
from v57_pie import read,write,sha
ROOT=old.ROOT;OUT=ROOT/'experiments/v80_hurdle_controls';A=OUT/'analysis'
ARMS=('hurdle_coupled','hurdle_independent','empirical')
data=old.data;contexts=old.contexts;key=old.key;checked_folder=old.checked_folder;bridge=old.bridge

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v80_hurdle_controls.py','v78_margin_variance.py','v29_methods.py','v29_pairwise_repair.py','v26_methods.py','v16_methods.py']
 write(A/'plan.json',dict(version=80,dependencies={n:sha(ROOT/n) for n in deps},arms=ARMS,hypothesis='Rank tests depend on cell CDFs and zero fractions, not only moments. Fit query NTC detection and positive-count shape; calibrate a relative NB detection response around empirical zeros, retaining binary control-cell dependence or testing independent detection. Empirical-scaled control transport is a comparator related to old V29, now at exact V66 realized margins. Specified without consulting V78/V79 DE summaries.',fit='Fit per-gene detection hazard against library offsets to observed average fraction nonzero (10 Newton moment iterations). Fit overdispersion of positive excess counts (X-1) conditional on positive cells/depth; shrink toward NTC alpha with npositive/(npositive+20). Use observed control bulk rate, alpha and median/mean depths; no query treatment input.',generation='Detection target uses logit empirical control detection plus logit NB-detection(mu_target,alpha)-logit NB-detection(mu_control,alpha). Cap detection<=mu_target to permit positive integer counts. Couple uniforms to actual control-template detection for coupled arm, independent uniforms otherwise. Given detection, sample1+NB(excess mean=mu_target/detection-1, fitted positive alpha). Condition at least one positive for a gene whose fixed parent total is positive. Exact V29 row/gene/CPM projection then enforces published V66 margins. Unsupported and current own genes receive no direct proposal; repair may affect carrier genes.',negative='Shuffle hazard, control alpha and positive alpha jointly among eligible genes within10 baseline-expression bins; couple against genuine template/control detection. Empirical arm permutes template gene columns within expression bins then rescales and projects. Norm matching of final cells not claimed.',constraints='Exact original row and gene counts => exact raw bulk profiles/PDS. Revert if CPM L1>.0005 or moved fraction>.35. Jackknife expression correction/DE means may change and are explicitly not protected.',selection='Same precision-focused V79 criterion, frozen before outcomes: development first2 unique cases/context; fidelity>=baseline+.01, Jaccard>=baseline+.001, reach>=baseline-.02, available-study NMAE<=1.01 and fidelity/Jaccard better than shuffled. Highest fidelity, tie Jaccard; no failed-guard fallback.',confirmation='Remaining20 cases, two seeds0/+620000, unchanged model. Both must pass critical criteria before production, exact profile/PDS and independent raw-p/rank checks. All public biology reused; not a hidden overall forecast or proof of true RNA kinetics.',D3_opened=False,automatic_upload=False))
 print('V80 empirical-zero hurdle protocol frozen',flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 old.check();return p

def fit_hurdle(x,q):
 x=sparse.csr_matrix(x,dtype=float);x.eliminate_zeros();n=x.shape[0];depth=np.asarray(x.sum(1)).ravel();meanlib=depth.mean();z=depth/meanlib;positive=x.copy();positive.data=np.ones_like(positive.data)
 npos=np.asarray(positive.sum(0)).ravel();fraction=npos/n;hazard=np.zeros(x.shape[1])
 for lo in range(0,x.shape[1],512):
  hi=min(lo+512,x.shape[1]);target=fraction[lo:hi];a=-np.log1p(-np.minimum(target,1-1e-8))
  for _ in range(10):
   zero=np.exp(-z[:,None]*a[None,:]);pred=1-zero.mean(0);derivative=(z[:,None]*zero).mean(0);a=np.clip(a+np.divide(target-pred,derivative,out=np.zeros_like(a),where=derivative>1e-12),0,100)
  hazard[lo:hi]=a
 excess=x.copy();excess.data-=1;excess.eliminate_zeros();total=np.asarray(excess.sum(0)).ravel();exposure=np.asarray(positive.T@depth).ravel();rate=np.divide(total,exposure,out=np.zeros_like(total),where=exposure>0);d2=np.asarray(positive.T@(depth*depth)).ravel();den=rate*rate*d2
 numerator=np.asarray(excess.power(2).sum(0)).ravel()-total-2*rate*np.asarray(excess.T@depth).ravel()+den
 alpha=np.divide(numerator,den,out=np.zeros_like(total),where=den>0);alpha=np.clip(alpha,0,20);weight=npos/(npos+20);alpha=weight*alpha+(1-weight)*q['alpha']
 return dict(hazard=hazard,alpha_positive=alpha,alpha=q['alpha'],p=q['p'],total=q['total'],mean_library=np.asarray(meanlib),genes=q['genes'],detection=fraction)

def prepare():
 check();from v25_methods import raw
 genes=data()['genes'];folder=A/'query_params';folder.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family'];path=folder/(name+'.npz')
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  x,meta,measured,source=raw(name,genes);ids=np.asarray(ctx['left']);assert (meta.gene_target.to_numpy(str)[ids]=='Non-Targeting').all();q=old.query(name);fit=fit_hurdle(x[ids],q)
  np.savez_compressed(path,**fit);write(path.with_suffix('.json'),dict(sha256=sha(path),source=name,controls=len(ids),source_record=source,script_sha256=sha(__file__)))
  print('V80 hurdle controls',name,len(ids),flush=True);del x;gc.collect()

def query(name):
 p=A/'query_params'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def gene_permutation(p,eligible):
 order=np.flatnonzero(eligible);order=order[np.argsort(p[order],kind='stable')];perm=np.arange(len(p))
 for block in np.array_split(order,10):
  if len(block)>1:perm[block]=np.roll(block,1)
 return perm

def nb_detection(mu,alpha):
 return -np.expm1(-np.log1p(alpha[None,:]*mu)/alpha[None,:])

def transform(parent,query,delta,arm,seed,shuffle=False,panel=None,template=None):
 parent=sparse.csr_matrix(parent,dtype=np.int32);original=parent.toarray();base=bridge.parameters(parent);libs=original.sum(1,dtype=np.int64);template=sparse.csr_matrix(template,dtype=np.int32).toarray()
 eligible=(query['total']>=80)&(base['total']>=20)
 if panel is not None:eligible&=~panel
 if not eligible.any():return parent.copy(),dict(CPM_L1=0.,fallback=False,changed_fraction=0.,bulk_exact=True,row_libraries_exact=True)
 perm=gene_permutation(query['p'],eligible) if shuffle else np.arange(len(eligible));rng=np.random.default_rng(seed);proposal=original.copy()
 if arm=='empirical':
  template=template[:,perm];template_mean=np.asarray(query['p'])[perm];scale=np.divide(base['p'],template_mean,out=np.ones_like(template_mean),where=template_mean>0)
  proposal=original.astype(float);proposal[:,eligible]=template[:,eligible]*scale[eligible]
 else:
  for lo in range(0,len(eligible),512):
   ix=np.arange(lo,min(lo+512,len(eligible)));ix=ix[eligible[ix]]
   if not len(ix):continue
   source=perm[ix];mu=libs[:,None]*base['p'][ix];mu0=libs[:,None]*query['p'][ix];z=libs[:,None]/float(query['mean_library'])
   baseline=-np.expm1(-z*query['hazard'][ix]);anchor=-np.expm1(-z*query['hazard'][source]);a=query['alpha'][source]
   before=nb_detection(mu0,a);after=nb_detection(mu,a)
   clip=lambda p:np.clip(p,1e-9,1-1e-9)
   d=expit(logit(clip(anchor))+logit(clip(after))-logit(clip(before)));d=np.minimum(d,mu);d=np.maximum(d,1e-12)
   u=rng.random(d.shape)
   if arm=='hurdle_coupled':u=np.where(template[:,ix]>0,u*baseline,baseline+u*(1-baseline))
   detected=u<d
   for j in np.flatnonzero(~detected.any(0)):
    weights=d[:,j];row=rng.choice(len(libs),p=weights/weights.sum());detected[row,j]=True
   positive_mean=np.maximum(mu/d-1,0);alpha=query['alpha_positive'][source];rate=positive_mean.copy();over=alpha>1e-10
   if over.any():rate[:,over]=rng.gamma(1/alpha[over],positive_mean[:,over]*alpha[over])
   proposed=(1+rng.poisson(rate))*detected;proposal[:,ix]=proposed
 count,diag=old.project(proposal,parent,seed+800080)
 np.testing.assert_array_equal(count.sum(0),parent.sum(0));np.testing.assert_array_equal(count.sum(1),parent.sum(1));moved=float(abs(count-parent).sum()/2/parent.sum());fallback=diag['mean_reference_L1']>.0005 or moved>.35
 if fallback:count=parent.copy()
 p0=np.asarray(parent.multiply((1/libs)[:,None]).mean(0)).ravel();p1=np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()
 return count,dict(diag,CPM_L1=float(abs(p1-p0).sum()),fallback=bool(fallback),changed_fraction=float(abs(count-parent).sum()/2/parent.sum()),bulk_exact=True,row_libraries_exact=True,reliable_genes=int(eligible.sum()))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='freeze':freeze()
  else:prepare()
