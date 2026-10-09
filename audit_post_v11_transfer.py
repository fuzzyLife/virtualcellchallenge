"""Held-family diagnostics of V11 uncertainty and target-specific transfer signal."""
import gc,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import read_cache
from v7_response_methods import FAMILIES,collapse,combine,adaptive_shrink
from v11_methods import consensus_gate
ROOT=build.ROOT;OUT=ROOT/'analysis/post_v11_research';SEEDS=[2026100901,2026100902]

def cosine(a,b,mask):
 a=np.where(mask,a,0);b=np.where(mask,b,0);den=np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
 return np.divide((a*b).sum(axis=1),den,out=np.full(len(a),np.nan),where=den>0)

def fold(data,family):
 train=data['row_families']!=family; y,var,valid,kd=collapse(data,train)
 gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 mean,variance,ok=combine((y,var,valid,kd),np.array([2.,.25,2.,1.,1.]));eb,_=adaptive_shrink(mean,variance,ok)
 gate,agreement,support=consensus_gate(y,var,valid,kd,[2.,.25,2.,1.,1.],mean,variance)
 pred=np.clip(eb+.2*gate*(mean-eb),-3,3)
 assert not np.any(valid[np.flatnonzero(FAMILIES==family)[0]])
 return pred,variance,ok,agreement,support

def main():
 OUT.mkdir(parents=True,exist_ok=True);data=read_cache();genes=data['genes'];targets=data['targets'];records=[];bins=[];coverage=[]
 _,_,valid,_=collapse(data,np.ones(len(data['truth']),bool));support=valid.sum(axis=0)
 for ti,target in enumerate(targets):
  assert support[ti,np.flatnonzero(genes==target)[0]]==0
  coverage.append({'target_gene':target,**{str(f)+'_usable_response_outputs':int(valid[fi,ti].sum()) for fi,f in enumerate(FAMILIES)},
   'outputs_any_source':int((support[ti]>0).sum()),'outputs_two_sources':int((support[ti]>=2).sum()),'outputs_no_source':int((support[ti]==0).sum()-1)})
 pd.DataFrame(coverage).to_csv(OUT/'source_target_coverage.csv',index=False);del valid,support
 for family in FAMILIES:
  pred,var,ok,agreement,support=fold(data,family)
  for ci in np.unique(data['row_contexts'][data['row_families']==family]):
   rows=np.flatnonzero((data['row_families']==family)&(data['row_contexts']==ci));ti=data['row_targets'][rows]
   truth=data['truth'][rows].astype(float);p=pred[ti].astype(float);mask=data['valid'][rows]&ok[ti]
   mask[np.arange(len(rows)),np.array([np.flatnonzero(genes==targets[t])[0] for t in ti])]=False
   de=mask&(data['q'][rows]<.05);err2=(p-truth)**2;scale=np.sqrt(var[ti]+data['errors'][rows]**2)
   z=np.divide(abs(p),scale,out=np.zeros_like(p),where=scale>0);wrong=np.sign(p)!=np.sign(truth)
   heuristic=norm.sf(z);cover=abs(p-truth)<=1.96*scale
   # Centering uses held truth solely as a descriptive evaluation transform, never fitting.
   diagnostic_rows=ok[ti].sum(axis=1)>=100
   common=ok[ti[diagnostic_rows]].all(axis=0)
   dmask=mask&common[None,:]&~np.isin(genes,targets)[None,:]
   dmask[~diagnostic_rows]=False
   n=dmask.sum(axis=0);pc=np.divide((p*dmask).sum(axis=0),n,out=np.zeros(len(genes)),where=n>0)
   tc=np.divide((truth*dmask).sum(axis=0),n,out=np.zeros(len(genes)),where=n>0)
   basecos=cosine(p,truth,dmask);rescos=cosine(p-pc,truth-tc,dmask)
   # A training-only pooled perturbation centroid baseline, excluding evaluated family.
   train=data['row_families']!=family;v=data['valid'][train];cnt=v.sum(axis=0)
   centroid=np.divide((data['truth'][train]*v).sum(axis=0),cnt,out=np.zeros(len(genes)),where=cnt>0)
   generic=cosine(np.broadcast_to(centroid,p.shape),truth,dmask)
   shuffled=[]
   for seed in SEEDS:
    active=np.flatnonzero(diagnostic_rows);perm=np.arange(len(rows));proposal=np.random.default_rng(seed+int(ci)).permutation(active)
    if np.any(proposal==active):proposal=np.roll(active,1+seed%(len(active)-1))
    perm[active]=proposal
    shuffled.append(cosine(pred[ti[perm]],truth,dmask))
   shuffled=np.asarray(shuffled);nsh=np.isfinite(shuffled).sum(axis=0)
   shuffle=np.divide(np.nansum(shuffled,axis=0),nsh,out=np.full(len(rows),np.nan),where=nsh>0)
   for j,ri in enumerate(rows):
    m=mask[j];dd=de[j];count=int(dd.sum());zero=float((truth[j,m]**2).mean()) if m.any() else np.nan
    records.append({'family':str(family),'context_index':int(ci),'target_gene':str(targets[ti[j]]),'evaluated_outputs':int(m.sum()),'diagnostic_common_outputs':int(dmask[j].sum()),
     'reference_de_normal_proxy':count,'lfc_mse':float(err2[j,m].mean()) if m.any() else np.nan,'zero_lfc_mse':zero,
     'lfc_mse_ratio':float(err2[j,m].mean()/zero) if zero>0 else np.nan,
     'de_sign_accuracy':float((~wrong[j,dd]).mean()) if count else np.nan,
     'heuristic_95_interval_coverage':float(cover[j,m].mean()) if m.any() else np.nan,'transfer_cosine':basecos[j],
     'residual_cosine_diagnostic':rescos[j],'shuffled_target_cosine':shuffle[j],'generic_train_centroid_cosine':generic[j]})
   for label,select in [('support1',support[ti]==1),('support2plus',support[ti]>=2),('agreement_low',agreement[ti]<.75),('agreement_high',agreement[ti]>=.75),
                        ('z0_1',z<1),('z1_2',(z>=1)&(z<2)),('z2_3',(z>=2)&(z<3)),('z3plus',z>=3)]:
    m=mask&select;dd=de&select
    bins.append({'family':str(family),'context_index':int(ci),'bin':label,'output_pairs':int(m.sum()),'de_pairs':int(dd.sum()),
     'heuristic_expected_sign_error_on_de':float(heuristic[dd].mean()) if dd.any() else np.nan,
     'observed_sign_error_on_de':float(wrong[dd].mean()) if dd.any() else np.nan,
     'heuristic_95_interval_coverage':float(cover[m].mean()) if m.any() else np.nan,
     'mean_squared_error':float(err2[m].mean()) if m.any() else np.nan})
   print('Audited transfer',family,int(ci),len(rows),'contrasts',flush=True)
  pd.DataFrame(records).to_csv(OUT/'transfer_target_diagnostics.csv',index=False);pd.DataFrame(bins).to_csv(OUT/'uncertainty_bins.csv',index=False)
  del pred,var,ok,agreement,support;gc.collect()
 df=pd.DataFrame(records);fields=[c for c in df if c not in ['family','context_index','target_gene']]
 summary=df.groupby(['family','context_index'])[fields].mean().groupby('family').mean();summary.to_csv(OUT/'transfer_family_summary.csv')
 build.write_json(OUT/'transfer_audit_manifest.json',{'script_sha256':build.sha256(Path(__file__)),
  'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
  'dependencies':{n:build.sha256(ROOT/n) for n in ['v7_response_methods.py','v11_methods.py','fit_v7_response.py']},
  'rows':len(df),'unsupported_contrasts':int((df.evaluated_outputs==0).sum()),'coverage_note':'Usable response masks include source abundance/QC filters; they are not complete assay gene axes. Own target excluded from no-source count.','seeds':SEEDS,'scope':'Cosine comparisons use the same common predictor-supported gene mask for matched and permuted targets in each context; unsupported contrasts abstain. Retrospective development diagnostic, previously inspected public source summaries. Entire evaluated family excluded. Normal-q source proxy, not count-level official metric. Ungated V11 estimator; no challenge-context gene-abundance gate. Heuristic intervals are not posterior credible intervals. Centered truth used only to calculate diagnostic metrics. No model selected or fitted for submission.'})
 print(summary[['lfc_mse_ratio','de_sign_accuracy','heuristic_95_interval_coverage','transfer_cosine','residual_cosine_diagnostic','shuffled_target_cosine','generic_train_centroid_cosine']].to_string(),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
