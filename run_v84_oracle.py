"""Diagnostic-only crossed mean/shape counts on disjoint donor treatment halves."""
import argparse,gc,json,hashlib
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v84_intervention_operator as m
from run_v84_donor import setup
from v57_pie import read,write,sha
from v16_methods import emit
from v6_methods import moments
from v29_methods import project
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from audit_v76_counts import spot_ranks

def run(chunk,chunks):
 protocol=read(m.A/'oracle_protocol.json');assessment=read(m.A/'donor_assessment_0.json');assert not assessment['passes']['program25'],'Conditional oracle not called on successful primary';x,meta,cohort,pred=setup();genes=pred['genes'];targets=pred['targets'];reference=x[cohort['right']];atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();selected=protocol['targets'];dest=m.A/'oracle';dest.mkdir(exist_ok=True)
 for j,target in enumerate(selected):
  if j%chunks!=chunk:continue
  case=next(c for c in cohort['cases'] if c['target']==target);view=next(v for v in case['views'] if v['offset']==0);fitted=x[view['actual'][:200]];actual=x[view['actual'][200:]];template=x[view['template'][:200]];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);fitlibs=np.asarray(fitted.sum(1)).ravel().astype(np.int64);own=list(genes).index(target);ti=list(targets).index(target);seed=view['seed']+8400000
  baseline,_=emit(template,np.arange(200),pred['v66'][ti],own,seed,'residual10',(libs,pred['cp'],pred['bulk']),atlas);shape_at_baseline,diag_shape=project(fitted.toarray()*np.divide(libs,fitlibs)[:,None],baseline,seed+1);old_at_measured,diag_mean=project(baseline.toarray()*np.divide(fitlibs,libs)[:,None],fitted,seed+2)
  np.testing.assert_array_equal(shape_at_baseline.sum(0),baseline.sum(0));np.testing.assert_array_equal(shape_at_baseline.sum(1),baseline.sum(1));np.testing.assert_array_equal(old_at_measured.sum(0),fitted.sum(0));np.testing.assert_array_equal(old_at_measured.sum(1),fitted.sum(1));truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
  for arm,count,diag in [('baseline',baseline,{}),('predicted_mean_fitted_shape',shape_at_baseline,diag_shape),('fitted_mean_old_shape',old_at_measured,diag_mean),('real_fit',fitted,{})]:
   p=dest/f'{target}__{arm}.json';pr,pj=profile_jackknife_chunked(count);dd=de(count,reference,genes);scores=direction_metrics(dd,truth,genes,own);spot=spot_ranks(count,reference,dd);axis=np.arange(len(genes))!=own
   arr=p.with_suffix('.npz');np.savez_compressed(arr,pred=pr,truth=tp,pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible']);write(p,dict(target=target,arm=arm,**scores,profile_MSE=float(np.mean((pr[axis]-tp[axis])**2)),projection_diagnostics=diag,CPM_contract_valid=diag.get('mean_reference_L1',0)<=.0005,manual_rank_coordinates_checked=spot,arrays_sha256=sha(arr),fit_cells=200,actual_cells=200,fit_and_actual_cells_disjoint=True,oracle_only=True,deployable=False));print('V84 oracle',target,arm,round(float(np.mean((pr[axis]-tp[axis])**2)),6),scores['jaccard'],flush=True)
  del baseline,shape_at_baseline,old_at_measured;gc.collect()

def assess():
 protocol=read(m.A/'oracle_protocol.json');records=[read(p) for p in (m.A/'oracle').glob('*.json')];assert len(records)==len(protocol['targets'])*4;summary={}
 for arm in sorted({r['arm'] for r in records}):
  rr=[r for r in records if r['arm']==arm];summary[arm]={}
  for k in ['profile_MSE','fidelity','jaccard','reach','nmae']:
   v=[r[k] for r in rr if r.get(k) is not None and np.isfinite(r[k])];summary[arm][k]=float(np.mean(v)) if v else None
  summary[arm]['CPM_contract_valid']=all(r['CPM_contract_valid'] for r in rr)
 for target in protocol['targets']:
  with np.load(m.A/'oracle'/f'{target}__baseline.npz') as b,np.load(m.A/'oracle'/f'{target}__predicted_mean_fitted_shape.npz') as q:np.testing.assert_array_equal(b['pred'],q['pred'])
  with np.load(m.A/'oracle'/f'{target}__real_fit.npz') as b,np.load(m.A/'oracle'/f'{target}__fitted_mean_old_shape.npz') as q:np.testing.assert_array_equal(b['pred'],q['pred'])
 write(m.A/'oracle_assessment.json',dict(summary=summary,records=records,raw_profile_identity_verified_within_each_mean_pair=True,oracle_only=True,query_treatment_labels_cannot_enter_inference=True,scope=protocol['limitations'],auditor='Per-population eight raw-rank coordinates; exact original gene/row margins and profile equality. This diagnostic is not a deployable model or final-test estimate.'));print(json.dumps(summary),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['run','assess']);p.add_argument('--chunk',type=int,default=0);p.add_argument('--chunks',type=int,default=1);a=p.parse_args()
 with threadpool_limits(limits=2):(run(a.chunk,a.chunks) if a.stage=='run' else assess())
