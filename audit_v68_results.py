"""Independently audit V68 cell selection, raw ranks and saved DE arithmetic."""
import gc,math
from pathlib import Path
import numpy as np
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
import v68_population_balancing as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from v62_conserved_response import data


def close(a,b,tol=1e-10):
    assert math.isclose(a,b,rel_tol=tol,abs_tol=tol),(a,b)


def independent_direction(z,genes,own):
    valid=z['truth_eligible']&(np.arange(len(genes))!=own)
    real=valid&(z['truth_q']<.05);adjudicable=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0)
    called=adjudicable&(z['pred_q']<.05);correct=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']))
    budget=max(np.count_nonzero(real),np.count_nonzero(called));union=real|(valid&(z['pred_q']<.05))
    fid=np.count_nonzero(correct)/budget if budget else None
    jac=np.count_nonzero(real&(z['pred_q']<.05))/np.count_nonzero(union) if union.any() else 1.
    pool=list(np.flatnonzero(real&adjudicable))
    ordered=sorted(pool,key=lambda i:(bool(z['pred_q'][i]>=.05),float(z['pred_q'][i]),float(z['pred_p'][i]),-abs(float(z['pred_lfc'][i])),str(genes[i])))
    reach=None
    if ordered:
        correct_so_far=0;best=0
        for j,i in enumerate(ordered,1):
            correct_so_far+=int(np.sign(z['pred_lfc'][i])==np.sign(z['truth_lfc'][i]))
            if correct_so_far/j>=.9:best=j
        reach=best/np.count_nonzero(real)
    elif real.any():reach=0.
    denom=float(abs(z['truth_lfc'][real]).sum())
    nmae=float(abs(z['pred_lfc'][real]-z['truth_lfc'][real]).sum())/denom if real.sum()>=10 and np.isfinite(denom) and denom>0 else None
    return dict(fidelity=fid,jaccard=jac,reach=reach,nmae=nmae)


def rank_p(a,b):
    n=len(a);k=len(b);combined=np.r_[a,b];r=stats.rankdata(combined,method='average');u=float(r[:n].sum()-n*(n+1)/2)
    _,ties=np.unique(combined,return_counts=True);total=n+k
    variance=n*k/12*((total+1)-np.sum(ties.astype(float)**3-ties)/(total*(total-1)))
    return float(2*stats.norm.sf(abs(u-n*k/2)/np.sqrt(variance))) if variance>0 else 1.


def main():
    m.check();d=data();genes=d['genes'];targets=d['targets'];receipts=[];sourcefiles=[];rank_checks=0
    for ctx in contexts():
        name=ctx['family'];folder,_=checked_folder(name);reference=sparse.load_npz(folder/'reference.npz')
        rf=m.A/'frames'/name;fr=read(rf/'receipt.json');assert sha(rf/'controls.npz')==fr['controls_sha256'];assert sha(rf/'cells.csv')==fr['cells_sha256']
        frame=sparse.load_npz(rf/'controls.npz');libs=frame.sum(1).A1;sourcefiles += [rf/'receipt.json',rf/'cells.csv']
        for path in sorted((m.A/'counts').glob(name+'__*.json')):
            if path.name.endswith('__panel.json'):continue
            r=read(path);sourcefiles.append(path)
            if r['arm']=='v63':
                prior=m.ROOT/r['reused_source'];assert sha(prior)==r['reused_receipt_sha256'];assert sha(prior.with_suffix('.npz'))==r['arrays_sha256'];continue
            arraypath=path.with_suffix('.npz');assert sha(arraypath)==r['arrays_sha256']
            with np.load(arraypath,allow_pickle=False) as f:z={k:f[k].copy() for k in f.files}
            case=ctx['cells'][r['case_index']];ti=list(targets).index(case['target']);own=list(genes).index(case['target']);seed=case['seed']+ti
            ids=z['selected_frame_rows'];draw=r['cells'];assert len(ids)==draw and len(np.unique(ids))==draw
            group,quota=m.strata(libs,draw);np.testing.assert_array_equal(np.bincount(group[ids],minlength=len(quota)),quota)
            edited=frame.copy();ownmask=edited.indices==own
            edited.data[ownmask]=np.random.default_rng(seed+680068).binomial(edited.data[ownmask],.1);edited.eliminate_zeros();c=edited[ids]
            difference=(frame[ids]-c).tocoo();assert (difference.col==own).all() and (difference.data>=0).all()
            np.testing.assert_array_equal(c.sum(1).A1,z['final_libraries']);np.testing.assert_array_equal(libs[ids],z['original_libraries'])
            pp=np.log1p(50000*c.sum(0).A1/float(c.sum()));np.testing.assert_allclose(pp,z['pred'],rtol=0,atol=2e-14)
            actual=sparse.load_npz(folder/f'{r["case_index"]}_actual.npz');truth_profile=np.log1p(50000*actual.sum(0).A1/float(actual.sum()))
            np.testing.assert_allclose(truth_profile,z['truth'],rtol=0,atol=2e-14)
            axis=np.arange(len(genes))!=own;close(float(np.sum((pp[axis]-truth_profile[axis])**2)),r['profile_square_error'])
            metrics=independent_direction(z,genes,own)
            for metric,value in metrics.items():
                if value is None:assert r[metric] is None
                else:close(value,r[metric])
            cl=c.sum(1).A1;rl=reference.sum(1).A1
            normalized=c.astype(float);normalized.data=np.divide(normalized.data,np.repeat(cl/1e6,np.diff(normalized.indptr)))
            normalized_ref=reference.astype(float);normalized_ref.data=np.divide(normalized_ref.data,np.repeat(rl/1e6,np.diff(normalized_ref.indptr)))
            lfc=np.log2((normalized.mean(0).A1+1e-9)/(normalized_ref.mean(0).A1+1e-9))
            np.testing.assert_allclose(lfc,z['pred_lfc'],atol=1e-12,rtol=1e-12)
            eligible=np.flatnonzero(z['pred_eligible']);selected=eligible[np.linspace(0,len(eligible)-1,min(16,len(eligible))).astype(int)]
            cl=c.sum(1).A1;rl=reference.sum(1).A1
            for g in selected:
                a=np.log1p(c[:,g].toarray().ravel()/(cl/1e6));b=np.log1p(reference[:,g].toarray().ravel()/(rl/1e6))
                close(rank_p(a,b),float(z['pred_p'][g]),tol=2e-10);rank_checks+=1
            # Independent BH reconstruction over all eligible genes.
            p=z['pred_p'][eligible];order=np.argsort(p,kind='stable');adjusted=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
            q=np.empty(len(p));q[order]=np.minimum(adjusted,1);np.testing.assert_allclose(q,z['pred_q'][eligible],atol=1e-13,rtol=1e-12)
            fitpath=m.A/'fits'/f'{name}__{r["case_index"]}.json';fitrec=read(fitpath);assert sha(fitpath)==r['fit_receipt_sha256'];assert sha(fitpath.with_suffix('.npz'))==fitrec['arrays_sha256'];sourcefiles.append(fitpath)
            with np.load(fitpath.with_suffix('.npz'),allow_pickle=False) as f:
                weights=f['weights'].copy();np.testing.assert_array_equal(group,f['group']);np.testing.assert_array_equal(quota,f['quota'])
            assert weights.min()>=0 and weights.max()<=1/draw+1e-12 and 1/(weights@weights)>=draw-1e-6
            for j,k in enumerate(quota):close(float(weights[group==j].sum()),k/draw,tol=1e-11)
            if r['arm']=='balanced_real':expected=m.sample(weights,group,quota,draw,seed+680000)
            elif r['arm']=='shuffled_weights':expected=m.sample(m.permute_within(weights,group),group,quota,draw,seed+680000)
            else:expected=m.sample(m.uniform(group,quota,draw),group,quota,draw,seed+680000)
            np.testing.assert_array_equal(ids,expected)
            receipts.append(dict(source=name,target=r['target'],arm=r['arm'],rows=draw,counts_reconstructed=True,non_own_counts_preserved=True,all_DE_metrics_recomputed=True))
        del frame;gc.collect();print('V68 audit completed',name,flush=True)
    assert len(receipts)==30
    for p in (m.A/'null').glob('*.json'):
        r=read(p);assert sha(p.with_suffix('.npz'))==r['arrays_sha256'];sourcefiles.append(p)
        with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:assert int(np.sum(z['eligible']&(z['q']<.05)))==r['calls']
    assert len(list((m.A/'null').glob('*.json')))==10
    result=read(m.A/'assessment.json');assert result['release_supported']==all(result['checks'].values())
    panels=[read(p) for p in (m.A/'counts').glob('*__panel.json')]
    for arm,summary in result['summary'].items():
        for metric,expected in summary.items():
            study_means=[]
            for study in ('H1','Orion','CD4'):
                vals=[p[metric] for p in panels if p['arm']==arm and p['study']==study and p.get(metric) is not None]
                study_means.append(sum(vals)/len(vals) if vals else None)
            if None in study_means:assert expected is None
            else:close(sum(study_means)/3,expected)
    complete_null=read(m.A/'primary_null_assessment.json');assert complete_null['passed']==all(r['passed'] for r in complete_null['records'])
    for path in (m.A/'primary_null').glob('*.json'):
        r=read(path);assert sha(path.with_suffix('.npz'))==r['arrays_sha256'];sourcefiles.append(path)
        with np.load(path.with_suffix('.npz'),allow_pickle=False) as z:
            assert len(np.unique(z['selected_frame_rows']))==400
            assert int(np.sum(z['eligible']&(z['q']<.05)))==r['calls']
    sourcefiles += [m.A/'primary_null_plan.json',m.A/'primary_null_assessment.json']
    sourcefiles += [m.A/'assessment.json',m.A/'plan.json',Path(__file__)]
    write(m.A/'results_audit.json',dict(passed=True,records=receipts,independent_raw_rank_coordinates=rank_checks,
          source_sha256={str(p.relative_to(m.ROOT)):sha(p) for p in sorted(set(sourcefiles))},
          scope='All30 new treatment populations reconstructed from original real-control rows and own-only edits; profiles, depth quotas, four DE metrics, all BH q-values and16 manual rank coordinates/population checked. Ten cached V63 references verified by hashes and frozen-run truth checks; no independent raw-p reach replay claimed for old caches. Null call totals reconciled; generator-null selected row profiles are not independently replayed here. Not official anchored accuracy.'))
    print('V68 independent count/rank audit passed:',len(receipts),'populations;',rank_checks,'manual rank coordinates.',flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()


