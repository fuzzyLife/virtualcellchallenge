"""V17 cross-fitted intervention-strength heterogeneity; V16 mean responses fixed."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
from scipy import sparse
import build_submission as b
from fit_v7_response import raw_inputs
from v16_methods import emit as baseline_emit
from v6_methods import moments
ROOT=b.ROOT;OUT=ROOT/'experiments/v17_learned_heterogeneity';A=OUT/'analysis'
RAW_FAMILIES=['H1','HCT116','HEK293T']
CANDIDATES=['v16','learned_half','learned_full','generic','shuffled']

def logcpm(x):
 depth=np.asarray(x.sum(1)).ravel();assert (depth>0).all()
 y=x.multiply((1e6/depth)[:,None]).tocsr();y.data=np.log1p(y.data);return y

def fit_source(family,genes,targets):
 x,meta,measured,rec=raw_inputs(family,genes);y=logcpm(x)
 panel=np.isin(genes,targets);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);batch_names=sorted(set(batches));bi={v:i for i,v in enumerate(batch_names)};batch_id=np.array([bi[v] for v in batches]);control={}
 for fold in [0,1]:
  means=[];variances=[];sizes=[]
  for sample in batch_names:
   ids=np.flatnonzero((labels=='Non-Targeting')&(batches==sample));ids=np.random.default_rng(17000+bi[sample]).permutation(ids)[fold::2]
   assert len(ids)>=8,'Insufficient matched controls'
   ct=y[ids];mu=np.asarray(ct.mean(0)).ravel();var=np.maximum(np.asarray(ct.power(2).mean(0)).ravel()-mu*mu,0);means.append(mu);variances.append(var*len(ids)/max(len(ids)-1,1));sizes.append(len(ids))
  control[fold]=(np.asarray(means),np.asarray(variances),np.asarray(sizes))
 records=[]
 for ti,target in enumerate(targets):
  ids=np.flatnonzero(labels==target)
  if len(ids)<80:continue
  rng=np.random.default_rng(np.random.SeedSequence([17001,RAW_FAMILIES.index(family),ti]));ids=rng.permutation(ids)
  for fold in [0,1]:
   tr=ids[fold::2];te=ids[1-fold::2];cm,cv,cn=control[fold];em,ev,en=control[1-fold]
   bw=np.bincount(batch_id[tr],minlength=len(batch_names))/len(tr);mu=bw@cm;var=bw@cv;delta=np.asarray(y[tr].mean(0)).ravel()-mu
   eligible=measured&~panel&(mu>.25)&(var>.01);signal=np.divide(abs(delta),np.sqrt(var+.1),out=np.zeros(len(genes)),where=eligible)
   selected=np.flatnonzero(eligible&(signal>.1));selected=selected[np.argsort(signal[selected])[-256:]]
   if len(selected)<10:continue
   w=delta[selected]/(var[selected]+.1);w/=np.linalg.norm(w)
   st=np.asarray(y[te][:,selected]@w).ravel()-em[batch_id[te]][:,selected]@w
   # Estimate matched within-batch CONTROL projection variance, preserving covariance.
   vc=0.;center_noise=0.
   for batch in np.unique(batch_id[te]):
    ctrl=np.flatnonzero((labels=='Non-Targeting')&(batch_id==batch));ctrl=np.random.default_rng(17000+int(batch)).permutation(ctrl)[1-fold::2]
    sc=np.asarray(y[ctrl][:,selected]@w).ravel();weight=float(np.mean(batch_id[te]==batch));vv=float(sc.var(ddof=1));vc+=weight*vv;center_noise+=weight*vv/len(ctrl)
   shift=float(st.mean());vt=float(st.var(ddof=1));noise=(vt+vc)/len(st);snr=shift*shift/max(noise,1e-12);cv2=max(vt-vc-center_noise,0)/max(shift*shift,1e-8);reliability=len(st)/(len(st)+100)*snr/(snr+10);estimate=float(min(cv2,.5)*reliability) if snr>=4 else 0.
   records.append({'family':family,'target':str(target),'target_index':ti,'sample':'matched_batch_centered_pool','fold':fold,'train_cells':len(tr),'test_cells':len(te),'control_train_cells':int(cn.sum()),'control_test_cells':int(en.sum()),'genes':len(selected),'response_shift':shift,'treated_projection_variance':vt,'control_projection_variance':vc,'control_centering_noise':center_noise,'snr':snr,'raw_cv2':cv2,'reliability':reliability,'cv2':estimate})
 pd.DataFrame(records).to_csv(A/f'{family}_variation_fits.csv',index=False)
 b.write_json(A/f'{family}_variation_provenance.json',{'family':family,'raw_input_record':rec,'methods_sha256':b.sha256(Path(__file__)),'rows':len(records),'training':'Only this source family, batch-matched centering; opposite treated/control cell halves fit response direction and projected dispersion; panel genes excluded. Subtract control-mean estimation noise.'})
 del x,y;gc.collect();print('V17 fitted variation source',family,len(records),flush=True)
 return records

def pool_variation(excluded,targets):
 frames=[pd.read_csv(A/f'{f}_variation_fits.csv') for f in RAW_FAMILIES if f not in excluded]
 fits=pd.concat(frames,ignore_index=True);means=fits.groupby(['family','target_index']).cv2.mean().groupby('target_index').mean()
 cv2=np.zeros(len(targets));cv2[means.index.to_numpy(int)]=means.to_numpy()
 assert np.isfinite(cv2).all() and ((cv2>=0)&(cv2<=.5)).all()
 generic=float(cv2.mean());shuffle=np.random.default_rng(17017).permutation(cv2)
 return {'cv2':cv2,'generic_cv2':np.full(len(targets),generic),'shuffled_cv2':shuffle,'support':np.bincount(fits.target_index.astype(int),minlength=len(targets))}

def strengths(model,name,ti):
 if name=='v16':return 0.
 if name=='generic':return float(model['generic_cv2'][ti])
 if name=='shuffled':return float(model['shuffled_cv2'][ti])
 return float(model['cv2'][ti])*(.5 if name=='learned_half' else 1.)

def emit(x,rows,effect,own,seed,cv2,prepared=None,atlas=None):
 if cv2<=0 or not np.any(effect):return baseline_emit(x,rows,effect,own,seed,'residual10',prepared,atlas)
 if not 0<=cv2<=.5:raise ValueError('Unbounded response strength variance')
 depth,p,bulk=prepared[:3] if prepared is not None else moments(x)[:3]
 rows=np.asarray(rows,int);template=x[rows];libs=np.asarray(depth)[rows].astype(np.int64)
 # Separate random stream leaves the Poisson stream unchanged at zero variance.
 rng=np.random.default_rng(np.random.SeedSequence([seed,17017]));z=rng.standard_normal(len(rows));z-=z.mean();z/=z.std()
 gain=1+np.sqrt(cv2)*np.clip(z,-1.4,2.8);gain=np.maximum(gain,.05);gain/=gain.mean()
 effect=np.clip(effect,-3,3);factor=np.exp2(effect);conditional=np.exp2(np.clip(gain[:,None]*effect[None,:],-3,3))
 if own is not None:factor[own]=.1;conditional[:,own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum()
 mean=p[None,:]*conditional;mean/=mean.sum(1,keepdims=True);expected=libs[:,None]*mean
 empirical=template.toarray().astype(float)*conditional;empirical*=np.divide(libs,empirical.sum(1))[:,None]
 initial=.25*expected+.65*np.random.default_rng(seed).poisson(expected)+.1*empirical;initial/=initial.sum(1,keepdims=True)
 if atlas is None:
  from evaluate_atlas_transfer import load_atlas_model
  atlas,_=load_atlas_model()
 count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
 np.testing.assert_array_equal(count.sum(1,dtype=np.int64),libs);assert np.issubdtype(count.dtype,np.integer) and (count>=0).all()
 rp=(count/libs[:,None]).mean(0);rb=count.sum(0,dtype=np.int64)/libs.sum()
 diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'response_gain_variance':float(gain.var()),'learned_cv2':float(cv2),'changed_transcript_fraction':float(abs(count-template.toarray()).sum()/2/libs.sum()),'global_moment_refit':True}
 if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Moment constraints failed')
 return sparse.csr_matrix(count),diag
