"""Independent arithmetic audit of the non-deployable crossed mean/shape diagnostic."""
import numpy as np
from v57_pie import read,write,sha
import v84_intervention_operator as m
from audit_v76_counts import close

def main():
 protocol=read(m.A/'oracle_protocol.json');cohort=read(m.ROOT/'experiments/v83_context_expansion/analysis/donor/cohort.json')
 with np.load(m.A/'donor_predictions.npz',allow_pickle=False) as z:genes=z['genes'].copy()
 records=[]
 for target in protocol['targets']:
  view=next(v for c in cohort['cases'] if c['target']==target for v in c['views'] if v['offset']==0)
  assert len(set(view['actual'][:200]))==len(set(view['actual'][200:]))==200
  assert not set(view['actual'][:200])&set(view['actual'][200:])
  own=list(genes).index(target);profiles={}
  for arm in ['baseline','predicted_mean_fitted_shape','fitted_mean_old_shape','real_fit']:
   path=m.A/'oracle'/f'{target}__{arm}.json';r=read(path);arr=path.with_suffix('.npz');assert sha(arr)==r['arrays_sha256']
   with np.load(arr,allow_pickle=False) as z:
    profiles[arm]=z['pred'].copy();axis=np.arange(len(genes))!=own
    close(float(np.mean((z['pred'][axis]-z['truth'][axis])**2)),r['profile_MSE'])
    for prefix in ['pred','truth']:
     ix=np.flatnonzero(z[prefix+'_eligible']);pv=z[prefix+'_p'][ix];order=np.argsort(pv,kind='stable');adjusted=np.minimum.accumulate((pv[order]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1];q=np.empty(len(ix));q[order]=np.minimum(adjusted,1);np.testing.assert_allclose(q,z[prefix+'_q'][ix],rtol=0,atol=1e-12)
    valid=z['truth_eligible']&axis;real=valid&(z['truth_q']<.05);judged=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0);called=judged&(z['pred_q']<.05);correct=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan
    union=real|(valid&(z['pred_q']<.05));jac=(real&(z['pred_q']<.05)).sum()/union.sum() if union.any() else 1.
    den=abs(z['truth_lfc'][real]).sum();nmae=abs(np.nan_to_num(z['pred_lfc'][real])-z['truth_lfc'][real]).sum()/den if real.sum()>=10 and np.isfinite(den) and den>0 else np.nan
    candidates=np.flatnonzero(real&judged);order=sorted(candidates,key=lambda i:(z['pred_q'][i]>=.05,z['pred_q'][i],z['pred_p'][i],-abs(z['pred_lfc'][i]),str(genes[i])))
    if order:
     purity=np.cumsum(np.sign(z['pred_lfc'][order])==np.sign(z['truth_lfc'][order]))/np.arange(1,len(order)+1);accepted=np.flatnonzero(purity>=.9);reach=(accepted[-1]+1)/int(real.sum()) if len(accepted) else 0.
    else:reach=0. if real.sum()>0 else np.nan
    for key,value in [('fidelity',fid),('jaccard',jac),('nmae',nmae),('reach',reach)]:close(value,r[key])
   assert r['manual_rank_coordinates_checked']==8 and r['oracle_only'] and not r['deployable'] and r['fit_and_actual_cells_disjoint']
   records.append(dict(path=path.relative_to(m.ROOT).as_posix(),sha256=sha(path),arrays_sha256=sha(arr),CPM_contract_valid=r['CPM_contract_valid']))
  np.testing.assert_array_equal(profiles['baseline'],profiles['predicted_mean_fitted_shape']);np.testing.assert_array_equal(profiles['real_fit'],profiles['fitted_mean_old_shape'])
 assert len(records)==24
 write(m.A/'oracle_audit.json',dict(passed=True,comparison_records=len(records),manual_rank_coordinates=192,independent_DE_BH_profile_replay=True,within_mean_pair_profiles_exact=True,fit_scoring_halves_disjoint=True,all_CPM_contracts_passed=all(r['CPM_contract_valid'] for r in records),records=records,auditor_sha256=sha(__file__),oracle_protocol_sha256=sha(m.A/'oracle_protocol.json'),non_deployable=True))
 print('V84 independent oracle arithmetic audit passed',len(records),flush=True)
if __name__=='__main__':main()
