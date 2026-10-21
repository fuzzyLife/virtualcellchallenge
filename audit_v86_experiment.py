"""Verify actual source likelihoods, exclusion, kernels and posterior calculations."""
import argparse,gc
from unittest.mock import patch
import numpy as np
from scipy import stats
from scipy.special import logsumexp
from threadpoolctl import threadpool_limits
import v86_joint_posterior as m
from v57_pie import read,write,sha

def models():
    m.check();assert sha(m.ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==read(m.A/'plan.json')['source_cache_sha256'];d=m.data();records=[]
    query_names={'H1':'H1','Orion':'HCT116','CD4':'CD4_Rest_D1'}
    for outer in ['H1','Orion','CD4']:
        b=m.bank(outer);expected={str(f) for f in m.FAMILIES if m.GROUP[str(f)]!=outer};assert set(b['source_names'])==expected
        for shuffled in [False,True]:
            net=m.model(outer,shuffled);path=m.A/'models'/(outer+('_shuffled' if shuffled else '')+'.npz');r=read(path.with_suffix('.json'));np.testing.assert_array_equal(net['similarity'],b['similarity']);np.testing.assert_array_equal(net['rank_features'],b['rank_features']);np.testing.assert_array_equal(net['axis'],b['axis']);y,v,ok=m.observations(b,shuffled);ids=net['sample_indices'];assert (ok[ids].sum(1)>=2).all();assert len(ids)==len(set(ids))==20000
            for weights in [net['weights'],net['shared_weights']]:assert (weights>=0).all();np.testing.assert_allclose(weights.sum(),1,atol=1e-12)
            ids=ids[:24];ll,_,_=m.likelihood(y[ids],v[ids],ok[ids],net['similarity'])
            for ci,(scale,length) in enumerate(m.COMPONENTS):
                if scale:assert np.linalg.eigvalsh(scale**2*m.kernel(net['similarity'],length)).min()>0
                for j,row in enumerate(ids):
                    use=np.flatnonzero(ok[row]);cov=np.diag(v[row,use].astype(float))
                    if scale:cov+=scale**2*m.kernel(net['similarity'],length)[np.ix_(use,use)]
                    manual=stats.multivariate_normal.logpdf(y[row,use],mean=np.zeros(len(use)),cov=cov);np.testing.assert_allclose(ll[j,ci],manual,rtol=1e-7,atol=1e-6)
            qp=m.A/'predictions'/(query_names[outer]+'.npz')
            with np.load(qp,allow_pickle=False) as z:query=z['query'].copy()
            axis=net['axis'];known=np.isfinite(query[axis]);feature=np.zeros(len(axis));rr=stats.rankdata(query[axis][known]);rr-=rr.mean();feature[known]=rr/np.linalg.norm(rr);sim=net['rank_features']@feature
            ll,means,var=m.likelihood(y[ids],v[ids],ok[ids],net['similarity'],sim)
            for ci,(scale,length) in enumerate(m.COMPONENTS):
                if not scale:continue
                for j,row in enumerate(ids[:8]):
                    use=np.flatnonzero(ok[row]);cov=np.diag(v[row,use].astype(float))+scale**2*m.kernel(net['similarity'],length)[np.ix_(use,use)];cross=scale**2*.8*(np.ones(len(use)) if length==999 else np.exp(-np.maximum(1-sim[use],0)/length));manual_mean=cross@np.linalg.solve(cov,y[row,use]);manual_var=scale**2-cross@np.linalg.solve(cov,cross);np.testing.assert_allclose(means[j,ci],manual_mean,rtol=1e-7,atol=1e-8);np.testing.assert_allclose(var[j,ci],manual_var,rtol=1e-7,atol=1e-8);assert var[j,ci]>=.2*scale**2-1e-8
            mu,pv=m.mixture(ll,means,var,net['weights']);p=np.exp(ll+np.log(np.maximum(net['weights'],1e-300))[None,:]-logsumexp(ll+np.log(np.maximum(net['weights'],1e-300))[None,:],axis=1)[:,None]);np.testing.assert_allclose(mu,np.sum(p*means,axis=1),atol=1e-12);np.testing.assert_allclose(pv,np.sum(p*(var+means**2),axis=1)-mu**2,atol=1e-12)
            records.append(dict(outer=outer,shuffled=shuffled,model_sha256=sha(path),fitting_coordinates=20000,likelihood_coordinates_checked=24,posterior_coordinates_checked=8,EM_iterations=len(r['em_history']),final_weight_change=r['em_history'][-1]['maximum_weight_change']))
            del y,v,ok;gc.collect()
        # Poison every recipient-study source response and profile before rebuilding.
        excluded=np.array([m.GROUP[str(f)]==outer for f in d['row_families']]);poison=dict(d)
        for key in ['truth','errors','valid','knockdown']:
            poison[key]=d[key].copy();poison[key][excluded]=True if key=='valid' else .01 if key=='knockdown' else 999.
        poison['profiles']=d['profiles'].copy();poison['profiles'][np.unique(d['row_contexts'][excluded])]=999999.
        m.bank.cache_clear()
        with patch.object(m,'data',return_value=poison):
            other=m.bank(outer)
            for key in ['y','variance','valid','source_profiles','rank_features','axis','similarity','source_indices','source_names']:np.testing.assert_array_equal(other[key],b[key])
        m.bank.cache_clear();del b,other,poison;gc.collect()
    prediction_records=[]
    for path in sorted((m.A/'predictions').glob('*.npz')):
        r=read(path.with_suffix('.json'));assert sha(path)==r['sha256']
        with np.load(path,allow_pickle=False) as z:
            base=z['v66'].copy();bulk=z['bulk'].copy();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);panel=np.isin(genes,d['targets']);control=np.log1p(50000*bulk);factor=np.exp2(base.astype(float));factor[np.arange(300),own]=.1;pp=factor*bulk;pp/=pp.sum(1,keepdims=True);delta=(np.log1p(50000*pp)-control)[:,~panel];norm=np.linalg.norm(delta,axis=1);u=np.divide(delta,norm[:,None],out=np.zeros_like(delta),where=norm[:,None]>0);s0=delta@u.T;floor=.8*np.maximum(np.diag(s0)[:,None]-s0,0);np.fill_diagonal(floor,-np.inf)
            for arm in m.ARMS:
                effect=z[arm];assert effect.shape==base.shape and np.isfinite(effect).all();np.testing.assert_allclose(effect[:,panel],base[:,panel],rtol=0,atol=5e-7);factor=np.exp2(effect.astype(float));factor[np.arange(300),own]=.1;pp=factor*bulk;pp/=pp.sum(1,keepdims=True);score=(np.log1p(50000*pp)-control)[:,~panel]@u.T;margin=np.diag(score)[:,None]-score;assert np.all(margin+1e-9*np.maximum(norm,1e-8)[:,None]>=floor)
        prediction_records.append(dict(path=path.relative_to(m.ROOT).as_posix(),sha256=sha(path)))
    assert len(prediction_records)==6
    write(m.A/'model_audit.json',dict(passed=True,models=records,predictions=prediction_records,whole_recipient_study_poison_has_no_effect=True,all_source_covariances_PSD=True,Gaussian_likelihoods_checked_against_scipy=True,conditional_posterior_checked_by_direct_solves=True,unseen_context_idiosyncratic_variance_retained=True,all_realized_identity_margins_replayed=True,auditor_sha256=sha(__file__)))
    print('V86 model, boundary and posterior audit passed',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['models','counts']);p.add_argument('--offset',type=int,default=0);a=p.parse_args()
    with threadpool_limits(limits=2):
        if a.stage=='models':models()
        else:
            from audit_v86_counts import main
            main(a.offset)
