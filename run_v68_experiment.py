"""Fixed small count/null test of V68's real-cell population-balancing hypothesis."""
import argparse,gc,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import v68_population_balancing as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from v25_methods import raw
from v62_conserved_response import data
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics

ARMS=('v63','balanced_real','shuffled_weights','uniform_real')


def unique_cases(ctx):
    seen=set();result=[]
    for i,c in enumerate(ctx['cells']):
        if c['target'] not in seen:result.append(i);seen.add(c['target'])
    return result[:2]


def prepare():
    m.check();d=data();genes=d['genes']
    for ctx in contexts():
        name=ctx['family'];folder=m.A/'frames'/name;folder.mkdir(parents=True,exist_ok=True);receipt=folder/'receipt.json'
        if receipt.exists():
            r=read(receipt);assert sha(folder/'controls.npz')==r['controls_sha256'] and sha(folder/'cells.csv')==r['cells_sha256'];continue
        x,meta,mask,source=raw(name,genes);left=np.asarray(ctx['left'],int);right=np.asarray(ctx['right'],int)
        assert len(set(left)&set(right))==0
        assert (meta.iloc[left].gene_target.astype(str)=='Non-Targeting').all()
        frame=x[left].astype(np.int64);assert np.array_equal(x[left].data,frame.data) and frame.data.min(initial=0)>=0
        libs=frame.sum(1).A1;assert (libs>0).all()
        field=next(f for f in ('cell_id','source_cell_id','cell') if f in meta)
        ids=meta.iloc[left][field].astype(str).to_numpy();assert len(set(ids))==len(ids)
        for i in unique_cases(ctx):
            assert not set(left)&set(ctx['cells'][i]['actual'])
        _,_,_,oldbulk,oldquery,_=reference_data(name)
        _,cp,bp,_=moments(frame)
        np.testing.assert_allclose(bp,oldbulk,rtol=1e-12,atol=1e-14)
        good=np.isfinite(oldquery);np.testing.assert_allclose(cp[good]*1e6,oldquery[good],rtol=1e-12,atol=1e-9)
        sparse.save_npz(folder/'controls.npz',frame)
        pd.DataFrame(dict(source_row=left,cell_id=ids,library=libs)).to_csv(folder/'cells.csv',index=False)
        write(receipt,dict(controls_sha256=sha(folder/'controls.npz'),cells_sha256=sha(folder/'cells.csv'),
              source=source,ntc_cells=len(left),reference_disjoint=True,treatments_disjoint=True,
              chosen_cases=unique_cases(ctx),plan_sha256=sha(m.A/'plan.json')))
        del x,frame;gc.collect();print('V68 cached real controls',name,len(left),flush=True)


def reference_data(name):
    folder,_=checked_folder(name);reference=sparse.load_npz(folder/'reference.npz')
    p=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz')
    assert sha(p)==read(p.with_suffix('.json'))['sha256']
    with np.load(p,allow_pickle=False) as z:base=z['v63'].copy();bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
    return folder,reference,base,bulk,query,measured


def score_counts(c,reference,truth,actual,genes,own):
    dd=de(c,reference,genes);pp,jp=profile_jackknife_chunked(c);tp,jt=profile_jackknife_chunked(actual)
    metrics=direction_metrics(dd,truth,genes,own);axis=np.arange(len(genes))!=own
    metrics['profile_square_error']=float(np.sum((pp[axis]-tp[axis])**2))
    return metrics,dict(pred=pp,jack=jp,truth=tp,truth_jack=jt,
          pred_p=dd['p'],pred_q=dd['q'],pred_lfc=dd['lfc'],pred_eligible=dd['eligible'],
          truth_p=truth['p'],truth_q=truth['q'],truth_lfc=truth['lfc'],truth_eligible=truth['eligible'])


def counts(only=None):
    m.check();d=data();genes=d['genes'];targets=d['targets'];dest=m.A/'counts';dest.mkdir(exist_ok=True)
    for ctx in contexts():
        name=ctx['family']
        if only is not None and name!=only:continue
        framepath=m.A/'frames'/name;fr=read(framepath/'receipt.json');assert sha(framepath/'controls.npz')==fr['controls_sha256']
        frame=sparse.load_npz(framepath/'controls.npz');original_lib=frame.sum(1).A1.astype(np.int64)
        folder,reference,base,bulk,query,measured=reference_data(name);cc,cj=profile_jackknife_chunked(reference)
        outputs={a:[] for a in ARMS}
        for i in unique_cases(ctx):
            cell=ctx['cells'][i];ti=list(targets).index(cell['target']);own=list(genes).index(cell['target']);seed=cell['seed']+ti
            actual=sparse.load_npz(folder/f'{i}_actual.npz');template=sparse.load_npz(folder/f'{i}_template.npz');draw=template.shape[0]
            truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
            prior=m.ROOT/'experiments/v66_assay_strength/analysis/counts'/f'{name}__{i}__0__v63.json';b=read(prior)
            assert sha(prior.with_suffix('.npz'))==b['arrays_sha256'] and b['cells']==draw and b['target']==cell['target']
            with np.load(prior.with_suffix('.npz'),allow_pickle=False) as z:
                old={a:z[a].copy() for a in z.files}
            np.testing.assert_array_equal(old['truth'],tp);np.testing.assert_array_equal(old['truth_q'],truth['q']);np.testing.assert_array_equal(old['truth_lfc'],truth['lfc'])
            axis=np.arange(len(genes))!=own
            baseline=dict(b,profile_square_error=float(np.sum((old['pred'][axis]-tp[axis])**2)),
                          reused_source=str(prior.relative_to(m.ROOT)),reused_receipt_sha256=sha(prior))
            write(dest/f'{name}__{i}__v63.json',baseline);outputs['v63'].append((old['pred'],old['jack'],tp,tj,baseline))
            edited=m.edit_own(frame,own,seed+680068);diff=(frame-edited).tocoo();assert (diff.col==own).all()
            group,quota=m.strata(original_lib,draw);factor=np.exp2(base[ti]);factor[own]=.1;target=bulk*factor;target/=target.sum()
            fitpath=m.A/'fits'/f'{name}__{i}.npz';fitpath.parent.mkdir(exist_ok=True)
            if fitpath.exists():
                fitrec=read(fitpath.with_suffix('.json'));assert sha(fitpath)==fitrec['arrays_sha256']
                with np.load(fitpath,allow_pickle=False) as z:p=z['weights'].copy();np.testing.assert_array_equal(z['group'],group);np.testing.assert_array_equal(z['quota'],quota)
            else:
                started=time.time();p,fitrec=m.fit(edited,target,measured&~np.isin(genes,targets),draw,group,quota)
                np.savez_compressed(fitpath,weights=p,group=group,quota=quota,target=target)
                fitrec=dict(fitrec,source=name,target=cell['target'],case_index=i,draw=draw,
                       arrays_sha256=sha(fitpath),frame_receipt_sha256=sha(framepath/'receipt.json'),seconds=time.time()-started)
                write(fitpath.with_suffix('.json'),fitrec)
            weights=dict(balanced_real=p,shuffled_weights=m.permute_within(p,group),uniform_real=m.uniform(group,quota,draw))
            for arm,w in weights.items():
                path=dest/f'{name}__{i}__{arm}.json';arrays=path.with_suffix('.npz')
                if not path.exists():
                    ids=m.sample(w,group,quota,draw,seed+680000);counts=edited[ids]
                    assert len(np.unique(ids))==draw and counts.data.min(initial=0)>=0
                    np.testing.assert_array_equal(counts.sum(1).A1,edited.sum(1).A1[ids])
                    result,values=score_counts(counts,reference,truth,actual,genes,own)
                    np.savez_compressed(arrays,**values,selected_frame_rows=ids,original_libraries=original_lib[ids],final_libraries=counts.sum(1).A1)
                    write(path,dict(family=name,study=key(name),case_index=i,target=cell['target'],arm=arm,cells=draw,
                           fit_receipt_sha256=sha(fitpath.with_suffix('.json')),arrays_sha256=sha(arrays),
                           unique_original_controls=True,non_own_counts_exactly_preserved=True,original_depth_quotas_exact=True,
                           effective_sample_size=float(1/(w@w)),**result))
                    print('V68 counts',name,cell['target'],arm,flush=True)
                r=read(path);assert sha(arrays)==r['arrays_sha256']
                with np.load(arrays,allow_pickle=False) as z:outputs[arm].append((z['pred'].copy(),z['jack'].copy(),z['truth'].copy(),z['truth_jack'].copy(),r))
            del edited;gc.collect()
        for arm,values in outputs.items():
            pp,jp,tt,jt=(np.asarray([v[k] for v in values]) for k in range(4));rows=[v[4] for v in values]
            panel=panel_metrics(pp,tt,cc,jp,jt,cj,genes,[v['target'] for v in rows],targets)
            for metric in ('fidelity','reach','jaccard','nmae','profile_square_error'):
                vals=[v[metric] for v in rows if v.get(metric) is not None and np.isfinite(v[metric])]
                panel[metric]=float(np.mean(vals)) if vals else None
            write(dest/f'{name}__{arm}__panel.json',dict(family=name,study=key(name),arm=arm,**panel))
        del frame;gc.collect()


def nulls():
    m.check();from v16_methods import emit
    from evaluate_atlas_transfer import load_atlas_model
    atlas,_=load_atlas_model();d=data();genes=d['genes'];dest=m.A/'null';dest.mkdir(exist_ok=True)
    for ci,ctx in enumerate(contexts()):
        name=ctx['family'];framepath=m.A/'frames'/name;fr=read(framepath/'receipt.json');assert sha(framepath/'controls.npz')==fr['controls_sha256']
        frame=sparse.load_npz(framepath/'controls.npz');folder,reference,base,bulk,query,measured=reference_data(name)
        group,quota=m.strata(frame.sum(1).A1,400);u=m.uniform(group,quota,400);seed=680800+ci
        for arm in ('v63','uniform_real'):
            output=dest/(name+'__'+arm+'.json');arrays=output.with_suffix('.npz')
            if output.exists():assert sha(arrays)==read(output)['arrays_sha256'];continue
            if arm=='uniform_real':
                ids=m.sample(u,group,quota,400,seed);counts=frame[ids]
            else:
                template=sparse.load_npz(folder/'null_template.npz');assert template.shape[0]==400
                libs=template.sum(1).A1.astype(np.int64)
                counts,_=emit(template,np.arange(400),np.zeros(len(genes)),None,seed,'residual10',(libs,np.nan_to_num(query)/1e6,bulk),atlas);ids=np.array([],int)
            dd=de(counts,reference,genes)
            np.savez_compressed(arrays,p=dd['p'],q=dd['q'],lfc=dd['lfc'],eligible=dd['eligible'],selected_frame_rows=ids)
            write(output,dict(family=name,arm=arm,cells=400,calls=int(np.sum(dd['eligible']&(dd['q']<.05))),arrays_sha256=sha(arrays)))
            print('V68 null',name,arm,read(output)['calls'],flush=True)
        del frame;gc.collect()


def assess():
    m.check();panels=[read(p) for p in (m.A/'counts').glob('*__panel.json')];assert len(panels)==20
    metrics=('pds','fidelity','reach','jaccard','nmae','profile_square_error');summary={};studies={}
    for arm in ARMS:
        studies[arm]={}
        for study in ('H1','Orion','CD4'):
            studies[arm][study]={}
            for metric in metrics:
                vals=[r[metric] for r in panels if r['arm']==arm and r['study']==study and r.get(metric) is not None and np.isfinite(r[metric])]
                studies[arm][study][metric]=float(np.mean(vals)) if vals else None
        summary[arm]={metric:float(np.mean([v[metric] for v in studies[arm].values()])) if all(v[metric] is not None for v in studies[arm].values()) else None for metric in metrics}
    b=summary['v63'];v=summary['balanced_real']
    def higher(k,other,loss=0):return v[k] is not None and other[k] is not None and v[k]>=other[k]-loss
    def improved(k,other):return v[k] is not None and other[k] is not None and (v[k]<other[k] if k=='nmae' else v[k]>other[k])
    ratio=float(np.mean([studies['balanced_real'][s]['profile_square_error']/studies['v63'][s]['profile_square_error'] for s in ('H1','Orion','CD4')]))
    checks=dict(profile=ratio<=.99,pds=higher('pds',b,.02),fidelity=higher('fidelity',b,.005),
                reach=higher('reach',b,.01),jaccard=higher('jaccard',b,.001),
                nmae=v['nmae'] is not None and b['nmae'] is not None and v['nmae']<=b['nmae']*1.01)
    de_gains=sum(improved(k,b) for k in ('fidelity','reach','jaccard','nmae'));checks['two_DE_gains']=de_gains>=2
    for arm in ('shuffled_weights','uniform_real'):
        other=summary[arm];checks['beats_'+arm+'_profile']=v['profile_square_error']<other['profile_square_error']
        gains=sum(improved(k,other) for k in ('fidelity','reach','jaccard','nmae'))
        checks['beats_'+arm+'_two_DE']=gains>=2
    null=[]
    for ctx in contexts():
        name=ctx['family'];a=read(m.A/'null'/(name+'__v63.json'));z=read(m.A/'null'/(name+'__uniform_real.json'))
        null.append(dict(source=name,v63_calls=a['calls'],natural_calls=z['calls'],passed=z['calls']<=max(5,a['calls'])))
    checks['nulls']=all(z['passed'] for z in null)
    write(m.A/'assessment.json',dict(summary=summary,studies=studies,checks=checks,nulls=null,
          paired_equal_study_profile_ratio=ratio,release_supported=all(checks.values()),
          unique_target_context_pairs=10,generated_comparisons=40,null_comparisons=10,
          no_private_accuracy_claim=True,D3_opened=False))
    print(json.dumps(dict(checks=checks,summary=summary,profile_ratio=ratio,nulls=null),indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','counts','nulls','assess']);p.add_argument('--source');a=p.parse_args()
    with threadpool_limits(limits=2):
        if a.stage=='counts':counts(a.source)
        else:globals()[a.stage]()

