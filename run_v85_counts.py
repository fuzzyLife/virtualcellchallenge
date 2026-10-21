"""V85 matched count diagnostic. Cohorts use metadata and source-supported predictions only."""
import argparse,gc,time,hashlib
import numpy as np
from threadpoolctl import threadpool_limits
from v57_pie import read,write,sha
import v85_guide_contrast as m
from v25_methods import raw
from run_v35_v36 import contexts
from v16_methods import emit
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics
from audit_v76_counts import spot_ranks

ARMS=('v66','contrast25','pooled25','target_permuted25')
NAMES=('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_D3_Rest')

def load(name):
    p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
    with np.load(p,allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
    if name=='CD4_D3_Rest':
        from run_v84_donor import setup
        x,meta,cohort,_=setup();left=cohort['left'];right=cohort['right'];receipt=read(m.ROOT/'experiments/v83_context_expansion/analysis/donor/extraction.json')
    else:
        x,meta,_,receipt=raw(name,m.data()['genes']);ctx=next(c for c in contexts() if c['family']==name);left=ctx['left'];right=ctx['right']
    assert not set(left)&set(right)
    assert (meta.iloc[left].gene_target=='Non-Targeting').all() and (meta.iloc[right].gene_target=='Non-Targeting').all()
    return x,meta,np.array(left),np.array(right),pred,receipt

def cohort():
    m.check();dest=m.A/'count_cohort.json';assert not dest.exists();cases=[];source_records=[]
    for name in NAMES:
        x,meta,left,right,pred,receipt=load(name);labels=meta.gene_target.to_numpy(str);lanes=meta['sample' if 'sample' in meta else 'lane_id'].to_numpy(str);left_by_lane={l:left[lanes[left]==l] for l in sorted(set(lanes[left]))};active=np.any(pred['contrast25']!=pred['v66'],axis=1);eligible=[]
        for ti,target in enumerate(m.data()['targets']):
            if not active[ti]:continue
            rng=np.random.default_rng(m.seed(name+str(target)+'cohort'));pool=[]
            for lane,controls in left_by_lane.items():
                rows=np.flatnonzero((labels==target)&(lanes==lane));pool.extend(rng.permutation(rows)[:len(controls)].tolist())
            if len(pool)<200:continue
            eligible.append((m.seed(name+str(target)+'case'),ti,str(target),np.array(pool)))
        for _,ti,target,pool in sorted(eligible,key=lambda r:r[0])[:2]:
            seed=m.seed(name+target+'population');rng=np.random.default_rng(seed);n=min(400,len(pool));actual=rng.choice(pool,n,replace=False);template=[]
            for lane in sorted(set(lanes[actual])):template.extend(rng.choice(left_by_lane[lane],int((lanes[actual]==lane).sum()),replace=False).tolist())
            rng.shuffle(template);assert len(set(template))==n and len(set(actual))==n
            cases.append(dict(source=name,target=target,target_index=ti,seed=seed,actual=actual.tolist(),template=template,n=n))
        source_records.append(dict(source=name,eligible_active_targets=len(eligible),reference=right.tolist(),fit_controls=left.tolist(),source_receipt=receipt,prediction_sha256=sha(m.A/'predictions'/(name+'.npz'))))
        print('V85 count cohort',name,'eligible',len(eligible),'retained',sum(c['source']==name for c in cases),flush=True);del x;gc.collect()
    write(dest,dict(cases=cases,sources=source_records,selected_before_count_scores=True,selection='At most2 smallest fixed source+target hashes among source-supported targets with metadata>=200 matched cells. No recipient response statistic used.',script_sha256=sha(__file__),plan_sha256=sha(m.A/'plan.json'),scope='Reused public studies and D3. D1/D3 use D2 new teacher explicitly as donor transfer. A small diagnostic, not independent biological confirmation.'))

def counts(chunk,chunks,offset):
    m.check();co=read(m.A/'count_cohort.json');atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();dest=m.A/'counts'/str(offset);dest.mkdir(parents=True,exist_ok=True);genes=m.data()['genes'];loaded=None
    for ci,case in enumerate(co['cases']):
        if ci%chunks!=chunk:continue
        name=case['source']
        if loaded!=name:
            x,meta,left,right,pred,receipt=load(name);reference=x[right];cp,cj=profile_jackknife_chunked(reference);loaded=name
        actual=x[case['actual']];template=x[case['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);target=case['target'];ti=case['target_index'];own=list(genes).index(target);truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
        for arm in ARMS:
            p=dest/f'{ci}__{arm}.json';arr=p.with_suffix('.npz')
            if p.exists():assert sha(arr)==read(p)['arrays_sha256'];continue
            start=time.time();count,diag=emit(template,np.arange(case['n']),pred[arm][ti],own,case['seed']+offset,'residual10',(libs,np.nan_to_num(pred['query'])/1e6,pred['bulk']),atlas)
            assert count.shape==(case['n'],len(genes)) and count.dtype.kind in 'iu' and (count.data>=0).all();np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
            pp,pj=profile_jackknife_chunked(count);dd=de(count,reference,genes);scores=direction_metrics(dd,truth,genes,own);axis=np.arange(len(genes))!=own;coordinates=spot_ranks(count,reference,dd)
            np.savez_compressed(arr,pred=pp,jack=pj,truth=tp,truth_jack=tj,control=cp,control_jack=cj,pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
            h=hashlib.sha256()
            for v in [count.data,count.indices,count.indptr]:h.update(np.ascontiguousarray(v).tobytes())
            write(p,dict(source=name,target=target,arm=arm,n=case['n'],case_index=ci,seed_offset=offset,profile_MSE=float(np.mean((pp[axis]-tp[axis])**2)),**scores,arrays_sha256=sha(arr),count_sha256=h.hexdigest(),original_NTC_libraries_exact=True,distinct_actual_and_template_cells=True,query_treatment_used_for_prediction=False,baseline_CD4_summary_excluded=name.startswith('CD4'),D2_teacher_allowed=True,manual_rank_coordinates_checked=coordinates,emitter=diag,seconds=time.time()-start));print('V85 counts',offset,name,target,arm,flush=True);del count;gc.collect()

def assess(offset):
    co=read(m.A/'count_cohort.json');genes=m.data()['genes'];targets=m.data()['targets'];dest=m.A/'counts'/str(offset);records=[];panels=[]
    for name in NAMES:
        cases=[(i,c) for i,c in enumerate(co['cases']) if c['source']==name]
        if not cases:continue
        for arm in ARMS:
            values=[];rows=[]
            for ci,case in cases:
                p=dest/f'{ci}__{arm}.json';r=read(p);assert sha(p.with_suffix('.npz'))==r['arrays_sha256'];rows.append(r);records.append(r)
                with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:values.append({k:z[k].copy() for k in ['pred','jack','truth','truth_jack','control','control_jack']})
            panel=panel_metrics(np.array([v['pred'] for v in values]),np.array([v['truth'] for v in values]),values[0]['control'],np.array([v['jack'] for v in values]),np.array([v['truth_jack'] for v in values]),values[0]['control_jack'],genes,[c['target'] for _,c in cases],targets) if len(cases)>=2 else dict(pds=None)
            for key in ['fidelity','reach','jaccard','nmae','profile_MSE']:
                finite=[r[key] for r in rows if r.get(key) is not None and np.isfinite(r[key])];panel[key]=float(np.mean(finite)) if finite else None
            panels.append(dict(source=name,study='Orion' if name in ['HCT116','HEK293T'] else 'CD4' if name.startswith('CD4') else name,arm=arm,cases=len(cases),**panel))
    summary={}
    for arm in ARMS:
        study={}
        for group in sorted(set(p['study'] for p in panels)):
            rows=[p for p in panels if p['study']==group and p['arm']==arm]
            study[group]={k:float(np.mean([r[k] for r in rows if r.get(k) is not None])) if any(r.get(k) is not None for r in rows) else None for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae']}
        summary[arm]={k:float(np.mean([v[k] for v in study.values()])) if study and all(v[k] is not None for v in study.values()) else None for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae']}
    b=summary['v66'];v=summary['contrast25'];checks={}
    for key,slack,smaller in [('profile_MSE',0,True),('pds',.005,False),('fidelity',.002,False),('reach',.005,False),('jaccard',.0005,False),('nmae',None,True)]:
        checks[key]=v[key] is not None and b[key] is not None and (v[key]<=b[key]*1.01 if key=='nmae' else v[key]<=b[key]+slack if smaller else v[key]>=b[key]-slack)
    gains=[v[k] is not None and b[k] is not None and (v[k]<b[k] if k in ['profile_MSE','nmae'] else v[k]>b[k]) for k in v];checks['three_gains']=sum(gains)>=3
    for control in ['pooled25','target_permuted25']:checks['beats_'+control]=v['profile_MSE'] is not None and summary[control]['profile_MSE'] is not None and v['profile_MSE']<summary[control]['profile_MSE']
    write(m.A/f'count_assessment_{offset}.json',dict(summary=summary,panels=panels,records=records,checks=checks,passes=all(checks.values()),seed_offset=offset,cohort_sha256=sha(m.A/'count_cohort.json'),scope='Small reused public diagnostic, unequal genuine n>=200, no anchor-scaled overall estimate. Equal study means; within CD4 and Orion average source contexts.'))
    print('V85 count assessment',offset,checks,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['cohort','counts','assess']);p.add_argument('--chunk',type=int,default=0);p.add_argument('--chunks',type=int,default=1);p.add_argument('--offset',type=int,default=0);a=p.parse_args()
    with threadpool_limits(limits=2):{'cohort':cohort,'counts':lambda:counts(a.chunk,a.chunks,a.offset),'assess':lambda:assess(a.offset)}[a.stage]()
