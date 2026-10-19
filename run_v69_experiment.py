"""V69 fixed source fitting, oracle diagnostic, transfer counts, selection and nulls."""
import argparse,gc,json,time
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v69_count_bridge as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from v62_conserved_response import data
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics


def unique(ctx):
    seen=set();out=[]
    for i,c in enumerate(ctx['cells']):
        if c['target'] not in seen:out.append(i);seen.add(c['target'])
    return out


def freeze():
    m.A.mkdir(parents=True,exist_ok=True);assert not (m.A/'plan.json').exists()
    deps=['v69_count_bridge.py','run_v69_experiment.py','test_v69_count_bridge.py',
          'v16_methods.py','v25_methods.py','v12_cpu_metrics.py','v6_methods.py','v8_methods.py',
          'audit_transfer_validation.py','run_v37_ensemble.py','run_v35_v36.py']
    inp=['experiments/v34_residual_flow/analysis/plan.json',
         'experiments/v66_assay_strength/analysis/reference_model.json',
         'experiments/v66_assay_strength/analysis/reference_model_arrays.npz']
    for c in contexts():
        inp.extend([f"experiments/v66_assay_strength/analysis/predictions/{c['family']}.npz",
                    f"experiments/v66_assay_strength/analysis/predictions/{c['family']}.json",
                    f"experiments/v68_population_balancing/analysis/frames/{c['family']}/receipt.json"])
    write(m.A/'plan.json',dict(version=69,dependencies={n:sha(m.ROOT/n) for n in deps},
        inputs={n:sha(m.ROOT/n) for n in inp},arms=m.ARMS,
        mechanism='Independent posterior Gamma-Poisson coupling with beta/gamma shape bridge, binomial removals and Poisson additions. Zero response exactly preserves counts. New cell totals may fluctuate; no row/gene integer-margin projection.',
        mean='Frozen whole-study-excluded V66 assay75 means; production exact frozen V66 effects. Own factor 0.1 applied once; target proportion closure includes own.',
        dispersion='Five raw source contexts, target minimum80 treated/80 matched NTC; per-gene depth-conditioned MOM alpha clipped .005..20. Split treated alternating source rows, same-sign dispersion shifts supported by molecule counts, shrink n/(n+200). Equal-study pooling plus one zero-change prior; transfer half clipped log alpha shift. All panel genes abstain for dispersion.',
        oracle='First2 unique targets/context; split cached actual cells alternating into fit and scoring halves. Five arms: old_pred, old_oracle, bridge_pred_oracle_shape, bridge_oracle, real_fit. All synthetic arms scored with same number of cells as scoring half. These are labelled-data diagnostics, not deployable predictors. Oracle old uses both measured fit CPM and bulk; bridge uses fit bulk/shape, so differences in realized means must be reported.',
        development={c['family']:unique(c)[:2] for c in contexts()},
        confirmation={c['family']:unique(c)[2:] for c in contexts()},
        selection='Select on development only among mean_bridge/shape_bridge/hybrid25: equal-study paired profile error ratio <=1.05, every study <=1.15, raw fidelity >=85% baseline, PDS no worse .03. Choose lowest paired profile error among eligible, tie by declared arm order. If none eligible, lock hybrid25 as explicitly experimental fallback, not a validated gain. Confirmation cannot change selected arm.',
        release_guards='Confirmation profile ratio <=.99, each study <=1.03, PDS no worse .005, fidelity no worse .002, reach no worse .005, Jaccard no worse .0005, NMAE <=1.01 baseline and >=3 raw metric gains. Beat matched shuffled shape control on profile for shape-based candidate. Complete zero-response and tiny-response null diagnostics reported separately. Failed guards require experimental label; user requested a new upload artifact, not automatic upload.',
        hybrid='Exactly one quarter of template row positions, selected by fixed seed permutation, replaced by full shape-bridge counts; remaining rows use unchanged old residual10 output. This mixture is not null-identity; no zero-response bypass is applied.',
        diagnostic_rule='Report whether bridge_oracle profile <=1.10 old_oracle and fidelity >=80% old_oracle; failure flags representation limitation. Still evaluate fixed deployable controls to meet the requested experimental release; do not claim oracle success or tune after seeing outcomes.',
        reserve='D3 remains unopened. All public biology reused; source study exclusion is not independent unseen-study confirmation.',
        automatic_upload=False))
    print('V69 protocol frozen',sha(m.A/'plan.json'),flush=True)


def frame(name):
    folder=m.ROOT/'experiments/v68_population_balancing/analysis/frames'/name
    r=read(folder/'receipt.json');assert sha(folder/'controls.npz')==r['controls_sha256']
    return sparse.load_npz(folder/'controls.npz')


def inputs(name):
    folder,_=checked_folder(name);path=m.ROOT/f'experiments/v66_assay_strength/analysis/predictions/{name}.npz'
    assert sha(path)==read(path.with_suffix('.json'))['sha256']
    with np.load(path,allow_pickle=False) as z:eff=z['assay75'].copy();query=z['query'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
    return folder,sparse.load_npz(folder/'reference.npz'),eff,query,bulk,measured


def score(path,count,reference,actual,genes,own,record):
    rawpath=path.with_suffix('.counts.npz');arrays=path.with_suffix('.npz')
    sparse.save_npz(rawpath,count);dd=de(count,reference,genes);truth=de(actual,reference,genes)
    pp,jp=profile_jackknife_chunked(count);tp,jt=profile_jackknife_chunked(actual)
    vals=direction_metrics(dd,truth,genes,own);axis=np.arange(len(genes))!=own
    vals['profile_square_error']=float(np.sum((pp[axis]-tp[axis])**2))
    vals['calls']=int(np.sum(dd['eligible']&(dd['q']<.05)))
    vals['mean_library']=float(count.sum(1).mean())
    np.savez_compressed(arrays,pred=pp,jack=jp,truth=tp,truth_jack=jt,p=dd['p'],q=dd['q'],lfc=dd['lfc'],eligible=dd['eligible'],truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])
    write(path,dict(record,**vals,arrays_sha256=sha(arrays),counts_sha256=sha(rawpath)))
    print('V69',record['stage'],record['family'],record['target'],record['arm'],'error',round(vals['profile_square_error'],4),'FID',vals['fidelity'],flush=True)


def old_emit(template,effect,own,seed,cp,bp):
    from v16_methods import emit
    from evaluate_atlas_transfer import load_atlas_model
    atlas,_=load_atlas_model();libs=template.sum(1).A1.astype(np.int64)
    return emit(template,np.arange(len(libs)),effect,own,seed,'residual10',(libs,cp,bp),atlas)[0]


def hybrid(old,new,seed):
    n=old.shape[0];ids=np.random.default_rng(seed+690069).permutation(n)[:n//4]
    out=old.tolil(copy=True);out[ids]=new[ids];return out.tocsr()


def run(stage,only=None):
    plan=m.check();d=data();genes=d['genes'];targets=d['targets'];dest=m.A/stage;dest.mkdir(exist_ok=True)
    selected=read(m.A/'selection.json')['selected'] if stage=='confirmation' else None
    for ctx in contexts():
        name=ctx['family']
        if only and only!=name:continue
        folder,reference,eff,q,bp,measured=inputs(name);controls=frame(name);params=m.parameters(controls)
        np.testing.assert_allclose(params['p'],bp,rtol=1e-12,atol=1e-14)
        delta,trace=m.pooled((key(name),));shuffled,_=m.pooled((key(name),),True)
        indices=unique(ctx)[:2] if stage in ('diagnostic','development') else unique(ctx)[2:]
        for i in indices:
            c=ctx['cells'][i];ti=list(targets).index(c['target']);own=list(genes).index(c['target']);seed=c['seed']+ti
            actual=sparse.load_npz(folder/f'{i}_actual.npz');template=sparse.load_npz(folder/f'{i}_template.npz')
            wanted=['old_pred','old_oracle','bridge_pred_oracle_shape','bridge_oracle','real_fit'] if stage=='diagnostic' else ['v66',*(m.ARMS if stage=='development' else (selected,)),'shuffled_shape']
            if all((dest/f'{name}__{i}__{arm}.json').exists() for arm in wanted):continue
            if stage=='diagnostic':
                fit=actual[::2];actual=actual[1::2];template=template[:actual.shape[0]]
                fitted=m.parameters(fit);od=np.log(fitted['alpha']/params['alpha'])
                oracle=np.log2(np.divide(fitted['p'],bp,out=np.ones_like(bp),where=(bp>0)&(fitted['p']>0)))
                oracle[(fitted['p']==0)&(bp>0)]=-3;oracle=np.clip(oracle,-3,3)
                populations={'old_pred':old_emit(template,eff[ti],own,seed,np.nan_to_num(q)/1e6,bp),
                             'old_oracle':old_emit(template,np.zeros(len(genes)),None,seed,moments(fit)[1],fitted['p']),
                             'bridge_pred_oracle_shape':m.bridge(template,eff[ti],own,seed,params,od,1.)[0],
                             'bridge_oracle':m.bridge(template,oracle,None,seed,params,od,1.)[0],
                             'real_fit':fit[:actual.shape[0]]}
            else:
                baseline=old_emit(template,eff[ti],own,seed,np.nan_to_num(q)/1e6,bp)
                shaped=m.bridge(template,eff[ti],own,seed+690000,params,delta[ti])[0]
                populations={'v66':baseline,'shape_bridge':shaped,'hybrid25':hybrid(baseline,shaped,seed),
                             'shuffled_shape':m.bridge(template,eff[ti],own,seed+690000,params,shuffled[ti])[0]}
                if stage=='development' or selected=='mean_bridge':populations['mean_bridge']=m.bridge(template,eff[ti],own,seed+690000,params)[0]
            for arm in wanted:
                path=dest/f'{name}__{i}__{arm}.json'
                if path.exists():continue
                count=populations[arm]
                score(path,count,reference,actual,genes,own,dict(stage=stage,family=name,study=key(name),case_index=i,target=c['target'],arm=arm,cells=count.shape[0],seed=seed,source_trace=trace,template_sha256=sha(folder/f'{i}_template.npz'),reference_sha256=sha(folder/'reference.npz'),actual_sha256=sha(folder/f'{i}_actual.npz'),plan_sha256=sha(m.A/'plan.json'),oracle=stage=='diagnostic'))
            del populations;gc.collect()
        panels(stage,ctx,wanted,reference,genes,targets)
        del controls;gc.collect()


def panels(stage,ctx,arms,reference,genes,targets):
    dest=m.A/stage;cp,cj=profile_jackknife_chunked(reference)
    for arm in arms:
        records=[read(p) for p in dest.glob(f"{ctx['family']}__*__{arm}.json")]
        if not records:continue
        records.sort(key=lambda r:r['case_index']);values=[]
        for r in records:
            path=dest/f"{ctx['family']}__{r['case_index']}__{arm}.npz";assert sha(path)==r['arrays_sha256']
            with np.load(path,allow_pickle=False) as z:values.append({k:z[k].copy() for k in ('pred','jack','truth','truth_jack')})
        metrics=panel_metrics(np.array([v['pred'] for v in values]),np.array([v['truth'] for v in values]),cp,np.array([v['jack'] for v in values]),np.array([v['truth_jack'] for v in values]),cj,genes,[r['target'] for r in records],targets)
        for metric in ('fidelity','reach','jaccard','nmae','profile_square_error'):
            vv=[r[metric] for r in records if r.get(metric) is not None and np.isfinite(r[metric])]
            metrics[metric]=float(np.mean(vv)) if vv else None
        write(dest/f"{ctx['family']}__{arm}__panel.json",dict(family=ctx['family'],study=key(ctx['family']),arm=arm,cases=len(records),**metrics))


def summary(stage):
    rows=[read(p) for p in (m.A/stage).glob('*__panel.json')];arms=sorted(set(r['arm'] for r in rows));ss={};av={}
    metrics=('pds','fidelity','reach','jaccard','nmae','profile_square_error')
    for arm in arms:
        ss[arm]={}
        for group in ('H1','Orion','CD4'):
            entries=[r for r in rows if r['arm']==arm and r['study']==group];assert entries,(stage,arm,group)
            ss[arm][group]={k:float(np.mean([r[k] for r in entries])) if all(r.get(k) is not None for r in entries) else None for k in metrics}
        av[arm]={k:float(np.mean([ss[arm][g][k] for g in ss[arm]])) if all(ss[arm][g][k] is not None for g in ss[arm]) else None for k in metrics}
    return dict(summary=av,studies=ss,scope='Raw reused-public metrics; not hosted anchor-scaled accuracy')


def assess(stage):
    m.check();s=summary(stage);av=s['summary'];ss=s['studies']
    if stage=='diagnostic':
        checks=dict(profile=av['bridge_oracle']['profile_square_error']<=1.1*av['old_oracle']['profile_square_error'],fidelity=av['bridge_oracle']['fidelity']>=.8*av['old_oracle']['fidelity'])
        write(m.A/'diagnostic_assessment.json',dict(s,checks=checks,representation_supported=all(checks.values()),not_a_deployable_predictor=True));print(json.dumps(dict(checks=checks,summary=av),indent=2),flush=True);return
    ratios={a:{g:ss[a][g]['profile_square_error']/ss['v66'][g]['profile_square_error'] for g in ss[a]} for a in av if a!='v66'}
    if stage=='development':
        eligible=[]
        for a in m.ARMS:
            if np.mean(list(ratios[a].values()))<=1.05 and max(ratios[a].values())<=1.15 and av[a]['fidelity']>=.85*av['v66']['fidelity'] and av[a]['pds']>=av['v66']['pds']-.03:eligible.append(a)
        selected=min(eligible,key=lambda a:np.mean(list(ratios[a].values()))) if eligible else 'hybrid25'
        path=m.A/'selection.json';assert not path.exists();write(path,dict(selected=selected,eligible=eligible,selection_stage='development_only',experimental_fallback=not bool(eligible),**s,ratios=ratios));print('V69 selected',selected,'eligible',eligible,flush=True)
    else:
        a=read(m.A/'selection.json')['selected'];v=av[a];b=av['v66']
        checks=dict(profile=np.mean(list(ratios[a].values()))<=.99,worst_study=max(ratios[a].values())<=1.03,pds=v['pds']>=b['pds']-.005,
             fidelity=v['fidelity']>=b['fidelity']-.002,reach=v['reach']>=b['reach']-.005,jaccard=v['jaccard']>=b['jaccard']-.0005,nmae=v['nmae']<=b['nmae']*1.01)
        gains=sum(v[k]>b[k] for k in ('pds','fidelity','reach','jaccard'))+int(v['nmae']<b['nmae'])+int(v['profile_square_error']<b['profile_square_error']);checks['three_gains']=gains>=3
        if a!='mean_bridge':checks['beats_shuffle']=v['profile_square_error']<av['shuffled_shape']['profile_square_error']
        write(m.A/'confirmation_assessment.json',dict(s,ratios=ratios,selected=a,checks=checks,release_supported=all(checks.values())));print(json.dumps(dict(checks=checks,summary=av,ratios=ratios),indent=2),flush=True)


def nulls(only=None):
    m.check();d=data();genes=d['genes'];zero=np.zeros(len(genes));dest=m.A/'null';dest.mkdir(exist_ok=True)
    for ci,ctx in enumerate(contexts()):
        name=ctx['family']
        if only and only!=name:continue
        folder,reference,_,q,bp,_=inputs(name);template=sparse.load_npz(folder/'null_template.npz');params=m.parameters(frame(name));seed=690800+ci
        old=old_emit(template,zero,None,seed,np.nan_to_num(q)/1e6,bp)
        exact,_=m.bridge(template,zero,None,seed,params)
        tiny=np.full(len(genes),1e-6);tiny[::2]*=-1
        near,_=m.bridge(template,tiny,None,seed,params)
        populations=dict(v66=old,bridge=exact,hybrid25=hybrid(old,exact,seed),near_zero_bridge=near)
        for arm,count in populations.items():
            path=dest/f'{name}__{arm}.json'
            if path.exists():continue
            dd=de(count,reference,genes);calls=int(np.sum(dd['eligible']&(dd['q']<.05)))
            write(path,dict(family=name,arm=arm,cells=count.shape[0],calls=calls,exact_count_identity=(count-template).nnz==0,changed_molecules=int(abs(count-template).sum()),count_sum=int(count.sum())))
            print('V69 null',name,arm,calls,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','diagnostic','development','confirmation','assess-diagnostic','select','assess-confirmation','null']);p.add_argument('--source');a=p.parse_args()
    if a.stage=='freeze':freeze()
    elif a.stage=='fit':
        m.check();d=data()
        for name in ([a.source] if a.source else m.SOURCES):m.fit_source(name,d['genes'],d['targets'])
    elif a.stage.startswith('assess-'):assess(a.stage.removeprefix('assess-'))
    elif a.stage=='select':assess('development')
    elif a.stage=='null':nulls(a.source)
    else:run(a.stage,a.source)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
