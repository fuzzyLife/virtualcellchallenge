"""V88 D3 count pilot with verified reuse of the identical V84 baseline."""
import argparse,gc,shutil,hashlib,time
import numpy as np
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import read,write,sha
from run_v84_donor import setup as old_setup
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked
from audit_v76_counts import spot_ranks
ARMS=('v66','exact125','shared125','exact125_shuffled')

def setup():
    m.check();x,meta,cohort,old=old_setup();path=m.A/'predictions/CD4_D3_Rest.npz';assert sha(path)==read(path.with_suffix('.json'))['sha256']
    with np.load(path,allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
    for key in ['v66','cp','bulk','genes','targets']:np.testing.assert_array_equal(pred[key],old[key])
    protocol=read(m.ROOT/'experiments/v84_intervention_operator/analysis/execution_protocol.json')
    for name in ['v16_methods.py','v12_cpu_metrics.py','v8_methods.py']:assert sha(m.ROOT/name)==protocol['code_sha256'][name]
    return x,meta,cohort,pred

def counts(offset,chunk,chunks):
    x,meta,cohort,pred=setup();genes=pred['genes'];targets=pred['targets'];reference=x[cohort['right']];cc,cj=profile_jackknife_chunked(reference);atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();dest=m.A/'counts'/str(offset);dest.mkdir(parents=True,exist_ok=True);reused=[]
    for ci,case in enumerate(cohort['cases']):
        if ci%chunks!=chunk:continue
        view=next(v for v in case['views'] if v['offset']==offset);target=case['target'];ti=list(targets).index(target);own=list(genes).index(target);actual=x[view['actual']];template=x[view['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
        source=m.ROOT/f'experiments/v84_intervention_operator/analysis/donor_counts/{offset}/{ci}__v66.json';record=read(source);arr=source.with_suffix('.npz');assert sha(arr)==record['arrays_sha256'];assert record['seed']==view['seed'] and record['target']==target
        with np.load(arr,allow_pickle=False) as z:
            np.testing.assert_array_equal(z['truth'],tp);np.testing.assert_array_equal(z['control'],cc);np.testing.assert_allclose(z['truth_lfc'],truth['lfc'],rtol=0,atol=0,equal_nan=True)
        target_json=dest/f'{ci}__v66.json';target_arr=target_json.with_suffix('.npz')
        for origin,output in [(source,target_json),(arr,target_arr)]:
            if output.exists():assert sha(output)==sha(origin)
            else:shutil.copyfile(origin,output)
        reused.append(dict(case_index=ci,source=source.relative_to(m.ROOT).as_posix(),source_sha256=sha(source),arrays_sha256=sha(arr)))
        for arm in ARMS[1:]:
            p=dest/f'{ci}__{arm}.json';arr=p.with_suffix('.npz')
            if p.exists():assert sha(arr)==read(p)['arrays_sha256'];continue
            start=time.time();count,diag=emit(template,np.arange(400),pred[arm][ti],own,view['seed'],'residual10',(libs,pred['cp'],pred['bulk']),atlas);assert count.shape==(400,len(genes)) and count.dtype.kind in 'iu' and (count.data>=0).all();np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
            pr,pj=profile_jackknife_chunked(count);dd=de(count,reference,genes);scores=direction_metrics(dd,truth,genes,own);coordinates=spot_ranks(count,reference,dd);axis=np.arange(len(genes))!=own;h=hashlib.sha256()
            for z in [count.data,count.indices,count.indptr]:h.update(np.ascontiguousarray(z).tobytes())
            np.savez_compressed(arr,pred=pr,jack=pj,truth=tp,truth_jack=tj,control=cc,control_jack=cj,pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
            write(p,dict(target=target,case_index=ci,arm=arm,seed_offset=offset,seed=view['seed'],prediction_cells=400,actual_cells=400,reference_cells=reference.shape[0],unique_prediction_control_cells=len(set(view['template'])),original_control_library_totals_exact=True,whole_CD4_response_family_excluded=True,query_treated_cells_used_for_inference=False,count_sha256=h.hexdigest(),profile_MSE=float(np.mean((pr[axis]-tp[axis])**2)),manual_rank_coordinates_checked=coordinates,**scores,emitter_diagnostics=diag,arrays_sha256=sha(arr),seconds=time.time()-start));print('V88 count',offset,target,arm,flush=True);del count;gc.collect()
    write(m.A/f'baseline_reuse_{offset}_{chunk}.json',dict(reused=reused,exact_V66_effect_CP_bulk_axis_and_cohort_match=True,original_emitter_and_metrics_code_match=True,truth_and_reference_profiles_exactly_replayed=True,new_baseline_cells_emitted=False,scope='Previously audited baseline counts reused; these are not new independent comparisons.'))

def assess(offset):
    _,_,cohort,pred=setup();dest=m.A/'counts'/str(offset);summary={};records=[]
    for arm in ARMS:
        rows=[];values=[]
        for ci,case in enumerate(cohort['cases']):
            p=dest/f'{ci}__{arm}.json';r=read(p);assert sha(p.with_suffix('.npz'))==r['arrays_sha256'];rows.append(r);records.append(r)
            with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:values.append({k:z[k].copy() for k in ['pred','jack','truth','truth_jack','control','control_jack']})
        panel=panel_metrics(np.array([v['pred'] for v in values]),np.array([v['truth'] for v in values]),values[0]['control'],np.array([v['jack'] for v in values]),np.array([v['truth_jack'] for v in values]),values[0]['control_jack'],pred['genes'],[r['target'] for r in rows],pred['targets'])
        for key in ['fidelity','reach','jaccard','nmae','profile_MSE']:
            vv=[r[key] for r in rows if r.get(key) is not None and np.isfinite(r[key])];panel[key]=float(np.mean(vv)) if vv else None
        summary[arm]=panel
    b=summary['v66'];v=summary['exact125'];checks={}
    for key,slack,smaller in [('profile_MSE',0,True),('pds',.005,False),('fidelity',.002,False),('reach',.005,False),('jaccard',.0005,False),('nmae',None,True)]:
        checks[key]=v[key] is not None and b[key] is not None and (v[key]<=1.01*b[key] if key=='nmae' else v[key]<=b[key]+slack if smaller else v[key]>=b[key]-slack)
    checks['three_gains']=sum(v[k] is not None and b[k] is not None and (v[k]<b[k] if k in ['profile_MSE','nmae'] else v[k]>b[k]) for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae'])>=3
    for control in ['shared125','exact125_shuffled']:checks['beats_'+control]=v['profile_MSE']<summary[control]['profile_MSE']
    write(m.A/f'count_assessment_{offset}.json',dict(summary=summary,checks=checks,passes=all(checks.values()),records=records,new_emitted_populations=42,reused_baseline_populations=14,seed_offset=offset,scope='Previously used D3 donor, whole CD4 excluded. Raw metrics, no hosted-score forecast or independent biology claim.'))
    print('V88 count assessment',offset,summary,checks,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['counts','assess']);p.add_argument('--offset',type=int,default=0);p.add_argument('--chunk',type=int,default=0);p.add_argument('--chunks',type=int,default=1);a=p.parse_args()
    with threadpool_limits(limits=2):(counts(a.offset,a.chunk,a.chunks) if a.stage=='counts' else assess(a.offset))
