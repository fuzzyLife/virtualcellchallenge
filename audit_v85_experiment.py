"""Independent source-statistic, exclusion, and count-metric audits for V85."""
import argparse,gc
from unittest.mock import patch
import numpy as np
from threadpoolctl import threadpool_limits
import v85_guide_contrast as m
from v57_pie import read,write,sha
from audit_v76_counts import close

def manual_estimate(y,var,own,base,nonpanel,measured):
    dose=y[:,own];dv=var[:,own];F=dose**2/np.maximum(dv,1e-12)
    valid=np.sign(dose[0])==np.sign(dose[1]) and min(abs(dose))>=.15 and min(F)>=10
    if not valid:return np.zeros(len(base)),np.zeros(len(base))
    a=-.9*y[1]/dose[0];b=-.9*y[0]/dose[1]
    va=.9**2*(np.sqrt(var[1])/abs(dose[0])+abs(y[1])*np.sqrt(dv[0])/dose[0]**2)**2
    vb=.9**2*(np.sqrt(var[0])/abs(dose[1])+abs(y[0])*np.sqrt(dv[1])/dose[1]**2)**2
    avg=(a+b)/2;se2=((np.sqrt(va)+np.sqrt(vb))/2)**2+((a-b)/2)**2
    uncertainty=se2/(np.log(2)*np.maximum(1+avg,.125))**2+.0025
    conf=.0625/(.0625+uncertainty)*(np.sign(a)==np.sign(b));effect=np.clip(np.log2(np.maximum(1+avg,.125)),-3,3)
    eligible=measured&nonpanel&(base>5);conf[~eligible]=0;effect[~eligible]=0
    if base[own]<=5:conf[:]=0
    return effect,conf

def sources():
    m.check();d=m.data();nonpanel=~np.isin(d['genes'],d['targets']);all_records=[]
    for condition in m.CONDITIONS:
        path=m.A/'linear_teachers'/(condition+'.npz');r=read(path.with_suffix('.json'));assert sha(path)==r['sha256'];sp=path.with_name(condition+'_statistics.npz');assert sha(sp)==r['statistics_sha256']
        with np.load(path,allow_pickle=False) as z:model={k:z[k].copy() for k in z.files}
        with np.load(sp,allow_pickle=False) as z:statistics={k:z[k].copy() for k in z.files}
        x,meta,measured,receipt=m.raw(condition,d['genes']);labels=meta.gene_target.to_numpy(str);lanes=meta['sample'].to_numpy(str);guides=meta.guide_id.to_numpy(str);depth=np.asarray(x.sum(1)).ravel()
        chosen={v['target'] for v in sorted(r['records'],key=lambda v:m.seed(v['target']))[:2]};replayed=[]
        def agg(rows):
            rows=np.array(rows);dense=x[rows].toarray().astype(float)*(1e6/depth[rows])[:,None]
            return dense.mean(0),dense.var(0,ddof=1)/len(rows)
        for prefix,key in [('', 'records'),('permuted_','permuted_records')]:
            for v in r[key]:
                ti=v['target_index'];s=str(ti);y=statistics[prefix+s+'_y'];var=statistics[prefix+s+'_variance'];base=statistics[prefix+s+'_base'];effect,conf=manual_estimate(y,var,v['own'],base,nonpanel,measured)
                np.testing.assert_allclose(effect,model[('permuted_effect' if prefix else 'effect')][ti],rtol=1e-6,atol=2e-7);np.testing.assert_allclose(conf,model[('permuted_confidence' if prefix else 'confidence')][ti],rtol=1e-6,atol=2e-8)
                flat=[]
                for lane in v['rows']:
                    ntc=lane['NTC_rows'];assert (labels[ntc]=='Non-Targeting').all() and (lanes[ntc]==lane['lane']).all()
                    for gi,halves in enumerate(lane['halves']):
                        for rows in halves:
                            assert (labels[rows]==v['target']).all() and (lanes[rows]==lane['lane']).all()
                            if not prefix:assert (guides[rows]==lane['guides'][gi]).all()
                            assert not set(rows)&set(ntc);flat.extend(rows)
                assert len(flat)==len(set(flat));assert len(flat)==2*sum(v['guide_half_cells'])
                if not prefix and v['target'] in chosen:
                    ds=[];vs=[];bs=[];cvs=[];weights=[]
                    for lane in v['rows']:
                        half=lane['halves'];group=[[agg(rows) for rows in h] for h in half];n1=len(half[0][0]);n2=len(half[1][0]);weights.append(n1*n2/(n1+n2));ds.append(np.array([group[0][h][0]-group[1][h][0] for h in [0,1]]));vs.append(np.array([group[0][h][1]+group[1][h][1] for h in [0,1]]));b,cv=agg(lane['NTC_rows']);bs.append(b);cvs.append(cv)
                    w=np.array(weights);w/=w.sum();ry=np.einsum('l,lhg->hg',w,np.array(ds));rv=np.einsum('l,lhg->hg',w*w,np.array(vs));bb=np.einsum('l,lg->g',w,np.array(bs));cv=np.einsum('l,lg->g',w*w,np.array(cvs));yy=np.divide(ry,bb[None,:],out=np.zeros_like(ry),where=bb[None,:]>0);vv=np.divide(rv,bb[None,:]**2,out=np.zeros_like(rv),where=bb[None,:]>0)+np.divide(ry**2*cv[None,:],bb[None,:]**4,out=np.zeros_like(ry),where=bb[None,:]>0)
                    np.testing.assert_allclose(yy,y,rtol=2e-7,atol=1e-8);np.testing.assert_allclose(vv,var,rtol=2e-7,atol=1e-8);np.testing.assert_allclose(bb,base,rtol=2e-7,atol=1e-8);replayed.append(v['target'])
        all_records.append(dict(condition=condition,source_sha256=sha(path),targets=len(r['records']),active_targets=r['active_targets'],permuted_active_targets=r['permuted_active_targets'],raw_statistic_replay_targets=replayed))
        del x,statistics,model;gc.collect()
    write(m.A/'source_audit.json',dict(passed=True,records=all_records,all_ratio_estimates_replayed=True,raw_count_normalization_replayed=True,all_source_cell_halves_disjoint=True,original_source_and_plan_preserved=True,auditor_sha256=sha(__file__)))
    print('V85 independent source audit passed',flush=True)

def predictions():
    m.check();d=m.data();axis=np.isin(d['genes'],d['targets']);records=[]
    for p in sorted((m.A/'predictions').glob('*.npz')):
        assert sha(p)==read(p.with_suffix('.json'))['sha256']
        with np.load(p,allow_pickle=False) as z:values={k:z[k].copy() for k in z.files}
        for arm in m.ARMS:
            np.testing.assert_array_equal(values[arm][:,axis],values['v66'][:,axis]);assert np.isfinite(values[arm]).all()
        with patch.object(m,'data',return_value={'genes':d['genes'],'targets':d['targets']}):
            again,_=m.predict(values['query'],values['v66'],values['bulk'],p.stem)
        for arm in m.ARMS:np.testing.assert_array_equal(again[arm],values[arm])
        with patch.object(m,'teachers',side_effect=AssertionError('Training donor attempted teacher access')):
            for recipient in m.CONDITIONS:
                none,_=m.predict(values['query'],values['v66'],values['bulk'],recipient)
                for arm in m.ARMS:np.testing.assert_array_equal(none[arm],values['v66'])
        records.append(dict(source=p.stem,sha256=sha(p)))
    assert len(records)==5
    write(m.A/'prediction_audit.json',dict(passed=True,records=records,recipient_response_cache_fields_removed_without_effect=True,training_donor_abstains_before_teacher_access=True,panel_and_own_effects_exactly_preserved=True,auditor_sha256=sha(__file__)))
    print('V85 prediction boundaries passed',flush=True)

def counts(offset):
    co=read(m.A/'count_cohort.json');folder=m.A/'counts'/str(offset);genes=m.data()['genes'];targets=m.data()['targets'];records=[];profiles={}
    for p in sorted(folder.glob('*.json')):
        r=read(p);arr=p.with_suffix('.npz');assert sha(arr)==r['arrays_sha256'];own=list(genes).index(r['target'])
        with np.load(arr,allow_pickle=False) as zz:z={k:zz[k].copy() for k in zz.files}
        for prefix in ['pred','truth']:
            ix=np.flatnonzero(z[prefix+'_eligible']);pv=z[prefix+'_p'][ix];order=np.argsort(pv,kind='stable');adjusted=np.minimum.accumulate((pv[order]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1];q=np.empty(len(ix));q[order]=np.minimum(adjusted,1);np.testing.assert_allclose(q,z[prefix+'_q'][ix],rtol=0,atol=1e-12)
        valid=z['truth_eligible']&(np.arange(len(genes))!=own);real=valid&(z['truth_q']<.05);judged=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0);called=judged&(z['pred_q']<.05);correct=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan;union=real|(valid&(z['pred_q']<.05));jac=(real&(z['pred_q']<.05)).sum()/union.sum() if union.any() else 1.
        den=abs(z['truth_lfc'][real]).sum();nmae=abs(np.nan_to_num(z['pred_lfc'][real])-z['truth_lfc'][real]).sum()/den if real.sum()>=10 and np.isfinite(den) and den>0 else np.nan
        ix=np.flatnonzero(real&judged);order=sorted(ix,key=lambda i:(z['pred_q'][i]>=.05,z['pred_q'][i],z['pred_p'][i],-abs(z['pred_lfc'][i]),str(genes[i])))
        if order:
            purity=np.cumsum(np.sign(z['pred_lfc'][order])==np.sign(z['truth_lfc'][order]))/np.arange(1,len(order)+1);ok=np.flatnonzero(purity>=.9);reach=(ok[-1]+1)/int(real.sum()) if len(ok) else 0.
        else:reach=0. if real.sum()>0 else np.nan
        for k,v in [('fidelity',fid),('jaccard',jac),('reach',reach),('nmae',nmae)]:close(v,r[k])
        axis=np.arange(len(genes))!=own;close(float(np.mean((z['pred'][axis]-z['truth'][axis])**2)),r['profile_MSE']);assert r['manual_rank_coordinates_checked']==8 and r['original_NTC_libraries_exact'] and not r['query_treatment_used_for_prediction'];case=co['cases'][r['case_index']];assert len(set(case['actual']))==len(set(case['template']))==r['n'];profiles.setdefault((r['source'],r['arm']),[]).append((r['case_index'],z['pred'],z['truth'],z['control']));records.append(dict(path=p.relative_to(m.ROOT).as_posix(),sha256=sha(p),arrays_sha256=sha(arr)))
    assessment=read(m.A/f'count_assessment_{offset}.json');axis=~np.isin(genes,targets)
    for panel in assessment['panels']:
        values=sorted(profiles[(panel['source'],panel['arm'])],key=lambda r:r[0])
        if len(values)<2:continue
        pp=(np.array([v[1] for v in values])-values[0][3])[:,axis];tt=(np.array([v[2] for v in values])-values[0][3])[:,axis];den=np.linalg.norm(pp,axis=1)[:,None]*np.linalg.norm(tt,axis=1)[None,:];sim=np.divide(pp@tt.T,den,out=np.zeros_like(den),where=den>0);ranks=[]
        for i in range(len(values)):
            other=np.delete(sim[i],i);ranks.append(np.mean((sim[i,i]>other)+.5*(sim[i,i]==other)))
        close(float(np.mean(ranks)),panel['pds'])
    assert len(records)==len(co['cases'])*4
    write(m.A/f'count_audit_{offset}.json',dict(passed=True,populations=len(records),manual_rank_coordinates=len(records)*8,independent_DE_BH_profile_PDS_replay=True,records=records,auditor_sha256=sha(__file__)))
    print('V85 independent count audit passed',len(records),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['sources','predictions','counts']);p.add_argument('--offset',type=int,default=0);a=p.parse_args()
    with threadpool_limits(limits=2):{'sources':sources,'predictions':predictions,'counts':lambda:counts(a.offset)}[a.stage]()
