"""Independent BH, direction metrics and target-rank replay for donor populations."""
import argparse
import numpy as np
from v57_pie import read,write,sha
import v89_projection_release as m
from audit_v76_counts import close

def main(offset):
 folder=m.A/'counts'/str(offset);records=[];gene_file=m.A/'predictions/CD4_D3_Rest.npz'
 with np.load(gene_file,allow_pickle=False) as z:genes=z['genes'].copy();targets=z['targets'].copy()
 for p in sorted(folder.glob('*.json')):
  r=read(p);arr=p.with_suffix('.npz');assert sha(arr)==r['arrays_sha256'];own=list(genes).index(r['target'])
  with np.load(arr,allow_pickle=False) as archive:
   z={k:archive[k].copy() for k in archive.files}
   valid=z['truth_eligible']&(np.arange(len(genes))!=own);real=valid&(z['truth_q']<.05);adjudicable=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0);called=adjudicable&(z['pred_q']<.05);correct=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan;union=real|(valid&(z['pred_q']<.05));jac=(real&(z['pred_q']<.05)).sum()/union.sum() if union.any() else 1.
   pool=np.flatnonzero(real&adjudicable);order=sorted(pool,key=lambda i:(z['pred_q'][i]>=.05,z['pred_q'][i],z['pred_p'][i],-abs(z['pred_lfc'][i]),str(genes[i])))
   if order:
    purity=np.cumsum(np.sign(z['pred_lfc'][order])==np.sign(z['truth_lfc'][order]))/np.arange(1,len(order)+1);ix=np.flatnonzero(purity>=.9);reach=(ix[-1]+1)/int(real.sum()) if len(ix) else 0.
   else:reach=0. if real.sum()>0 else np.nan
   den=abs(z['truth_lfc'][real]).sum();nmae=abs(np.nan_to_num(z['pred_lfc'][real])-z['truth_lfc'][real]).sum()/den if real.sum()>=10 and np.isfinite(den) and den>0 else np.nan
   for key,value in [('fidelity',fid),('jaccard',jac),('reach',reach),('nmae',nmae)]:close(value,r[key])
   axis=np.arange(len(genes))!=own;close(float(np.mean((z['pred'][axis]-z['truth'][axis])**2)),r['profile_MSE'])
   for prefix in ['pred','truth']:
    indices=np.flatnonzero(z[prefix+'_eligible']);pv=z[prefix+'_p'][indices];order=np.argsort(pv,kind='stable');manual=np.minimum.accumulate((pv[order]*len(indices)/np.arange(1,len(indices)+1))[::-1])[::-1];q=np.empty(len(indices));q[order]=np.minimum(manual,1);np.testing.assert_allclose(q,z[prefix+'_q'][indices],rtol=0,atol=1e-12)
  assert r['prediction_cells']==r['actual_cells']==r['unique_prediction_control_cells']==400;assert r['original_control_library_totals_exact'] and r['whole_CD4_response_family_excluded'] and not r['query_treated_cells_used_for_inference'];assert r['manual_rank_coordinates_checked']==8;records.append(dict(path=p.relative_to(m.ROOT).as_posix(),sha256=sha(p),arrays_sha256=sha(arr)))
 assessment=read(m.A/f'count_assessment_{offset}.json');panels=[];cohort=read(m.ROOT/'experiments/v83_context_expansion/analysis/donor/cohort.json')
 for arm,summary in assessment['summary'].items():
  values=[]
  for ci,case in enumerate(cohort['cases']):
   with np.load(folder/f'{ci}__{arm}.npz',allow_pickle=False) as z:values.append((z['pred'].copy(),z['truth'].copy(),z['control'].copy()))
  mask=~np.isin(genes,targets);pp=(np.array([v[0] for v in values])-values[0][2])[:,mask];tt=(np.array([v[1] for v in values])-values[0][2])[:,mask];den=np.linalg.norm(pp,axis=1)[:,None]*np.linalg.norm(tt,axis=1)[None,:];sim=np.divide(pp@tt.T,den,out=np.zeros_like(den),where=den>0);ranks=[]
  for i in range(len(values)):
   other=np.delete(sim[i],i);ranks.append(np.mean((sim[i,i]>other)+.5*(sim[i,i]==other)))
  value=float(np.mean(ranks));close(value,summary['pds']);panels.append(dict(arm=arm,pds=value))
 expected=len(cohort['cases'])*len(assessment['summary']);assert len(records)==expected
 write(m.A/f'count_audit_{offset}.json',dict(passed=True,comparison_records=len(records),manual_tied_rank_coordinates=len(records)*8,independent_DE_BH_profile_and_PDS_replay=True,records=records,panels=panels,seed0_known_at_selection=True,second_seed_not_used_for_postfreeze_hyperparameter_selection=True,new_populations_emitted=42,reused_baseline_populations=14,new_tied_rank_coordinates=336,reused_tied_rank_coordinates=112,auditor_sha256=sha(__file__),scope='Numerical replay of42 new and14 verified reused baseline populations. D3 was used in V84; it is not independent or blind biology. The112 baseline raw-rank coordinates were previously checked;336 new coordinates were checked during V89 emission.'))
 print('V89 independent count metric audit passed',offset,len(records),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--offset',type=int,default=0);main(p.parse_args().offset)
