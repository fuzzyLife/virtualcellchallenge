"""Broad, target-split response test for the frozen V67 intervention network."""
import argparse
import gc
import hashlib
import json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v67_intervention_network as m
from v57_pie import read, write, sha
from run_v35_v36 import contexts, key
from run_v37_ensemble import checked_folder
from v25_methods import raw
from audit_transfer_validation import profile, pds
from v6_methods import moments

ARMS=('v63','v66')+m.ARMS


def split(target):
    return 'development' if int(hashlib.sha256(('v67-target:'+target).encode()).hexdigest()[:8],16)%2==0 else 'confirmation'


def means():
    m.check();d=m.old.data();genes=d['genes'];targets=d['targets'];all_records=[];panels=[]
    dest=m.A/'broad';dest.mkdir(exist_ok=True)
    for ctx in contexts():
        name=ctx['family'];output=dest/(name+'.json')
        if output.exists():
            r=read(output);assert r['plan_sha256']==sha(m.A/'plan.json')
            assert sha(m.A/'predictions'/(name+'.npz'))==r['predictions_sha256']
            all_records+=r['records'];panels+=r['panels'];continue
        path=m.A/'predictions'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256']
        with np.load(path,allow_pickle=False) as z:pred={a:z[a].copy() for a in ARMS};bulk=z['bulk'].copy();query=z['query'].copy();measured=z['measured'].copy()
        folder,_=checked_folder(name);reference=sparse.load_npz(folder/'reference.npz');ref=profile(reference);ref_cpm=moments(reference)[1]*1e6
        x,meta,mask,source=raw(name,genes);labels=meta.gene_target.to_numpy(str)
        chosen=[(ti,str(t),np.flatnonzero(labels==t)) for ti,t in enumerate(targets) if np.count_nonzero(labels==t)>=30]
        assert len(chosen)>0
        tp=[];cpm=[]
        for ti,t,ids in chosen:
            tp.append(profile(x[ids]));cpm.append(moments(x[ids])[1]*1e6)
        tp=np.asarray(tp);cpm=np.asarray(cpm);records=[];ps=[]
        groups={'all':np.arange(len(chosen)), 'development':np.array([i for i,(_,t,_) in enumerate(chosen) if split(t)=='development']),
                'confirmation':np.array([i for i,(_,t,_) in enumerate(chosen) if split(t)=='confirmation']),
                'canonical_small':np.array([i for i,(_,t,_) in enumerate(chosen) if t in {c['target'] for c in ctx['cells']}])}
        scoring=measured&mask; pds_axis=scoring&~np.isin(genes,targets)
        for arm,effects in pred.items():
            pp=[]
            for j,(ti,t,ids) in enumerate(chosen):
                own=int(np.flatnonzero(genes==t)[0]);factor=np.exp2(effects[ti]);factor[own]=.1
                probability=bulk*factor;probability/=probability.sum();pr=np.log1p(50000*probability);pp.append(pr)
                axis=scoring.copy();axis[own]=False
                valid=axis&(ref_cpm>5)
                cp=np.nan_to_num(query)*factor;cp*=1e6/cp.sum()
                lfc=np.clip(np.log2((cp+.1)/(ref_cpm+.1)),-3,3)
                truth=np.clip(np.log2((cpm[j]+.1)/(ref_cpm+.1)),-3,3)
                records.append(dict(source=name,study=key(name),target=t,split=split(t),cells=len(ids),arm=arm,
                     canonical_small=t in {c['target'] for c in ctx['cells']},
                     treatment_rows_sha256=hashlib.sha256(ids.astype(np.int64).tobytes()).hexdigest(),
                     profile_MSE=float(np.mean((pr[axis]-tp[j,axis])**2)),
                     lfc_MSE=float(np.mean((lfc[valid]-truth[valid])**2)),
                     lfc_MAE=float(np.mean(abs(lfc[valid]-truth[valid])))))
            pp=np.asarray(pp)
            for group,idx in groups.items():
                ps.append(dict(source=name,study=key(name),arm=arm,group=group,targets=len(idx),
                     pds=pds(pp[idx][:,pds_axis],tp[idx][:,pds_axis],ref[pds_axis]) if len(idx)>1 else None))
        r=dict(plan_sha256=sha(m.A/'plan.json'),predictions_sha256=sha(path),
               raw_source=source,targets=len(chosen),target_names=[t for _,t,_ in chosen],records=records,panels=ps)
        write(output,r);all_records+=records;panels+=ps
        print('V67 broad response',name,len(chosen),flush=True);del x, pred;gc.collect()
    summaries={}
    for arm in m.ARMS:
        summaries[arm]={}
        for baseline in ('v63','v66'):
            summaries[arm][baseline]={}
            for group in ('all','development','confirmation','canonical_small'):
                studies={}
                for study in ('H1','Orion','CD4'):
                    src=sorted({r['source'] for r in all_records if r['study']==study})
                    ratios={k:[] for k in ('profile_MSE','lfc_MSE','lfc_MAE')};pchanges=[]
                    for name in src:
                        def sel(r,a):
                            return r['source']==name and r['arm']==a and (group=='all' or (r['canonical_small'] if group=='canonical_small' else r['split']==group))
                        for k in ratios:
                            aa=[r[k] for r in all_records if sel(r,arm)];bb=[r[k] for r in all_records if sel(r,baseline)]
                            assert len(aa)==len(bb) and aa
                            ratios[k].append(float(np.mean(aa)/np.mean(bb)))
                        pa=next(r['pds'] for r in panels if r['source']==name and r['arm']==arm and r['group']==group)
                        pb=next(r['pds'] for r in panels if r['source']==name and r['arm']==baseline and r['group']==group)
                        if pa is not None and pb is not None:pchanges.append(pa-pb)
                    studies[study]={k+'_ratio':float(np.mean(v)) for k,v in ratios.items()}
                    studies[study]['pds_change']=float(np.mean(pchanges))
                summaries[arm][baseline][group]=dict(studies=studies,**{k:float(np.mean([v[k] for v in studies.values()])) for k in next(iter(studies.values()))})
    checks={}
    for baseline in ('v63','v66'):
        s=summaries['network'][baseline]
        checks['profile_vs_'+baseline]=s['all']['profile_MSE_ratio']<=.99
        checks['worst_study_vs_'+baseline]=max(v['profile_MSE_ratio'] for v in s['all']['studies'].values())<=1.02
        checks['pds_vs_'+baseline]=s['all']['pds_change']>=-.005
        for group in ('development','confirmation'):
            checks[group+'_profile_vs_'+baseline]=s[group]['profile_MSE_ratio']<1
    for control in ('network_shuffled','gain_only'):
        checks['beats_'+control]=summaries['network']['v63']['all']['profile_MSE_ratio']<summaries[control]['v63']['all']['profile_MSE_ratio']
    result=dict(plan_sha256=sha(m.A/'plan.json'),summary=summaries,checks=checks,
                primary='network',selected='network' if all(checks.values()) else None,
                release_supported=all(checks.values()),unique_target_context_pairs=sum(read(dest/(c['family']+'.json'))['targets'] for c in contexts()),
                reused_public_biology=True,D3_opened=False,private_score=None,
                scope='Intended bulk means before count generation. Whole-study exclusion. Fixed target partitions reuse prior public biology and are not independent blind validation. Local diagnostics are not official anchored scores.')
    write(m.A/'mean_results.json',result)
    print(json.dumps(dict(checks=checks,summary={a:summaries[a]['v63']['all'] for a in m.ARMS}),indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['means']);args=p.parse_args()
    with threadpool_limits(limits=2):globals()[args.stage]()
