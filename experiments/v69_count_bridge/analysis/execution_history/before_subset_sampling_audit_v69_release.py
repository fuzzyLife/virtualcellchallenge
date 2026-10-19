"""Independent checks of V69 saved counts, DE arithmetic and deployed payload."""
from pathlib import Path
import argparse,gc,json,tempfile,subprocess
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
import v69_count_bridge as m
from v57_pie import read,write,sha
from build_v69_submission import load,SEED
from run_v35_v36 import contexts


def normalized(x):
    x=sparse.csr_matrix(x,dtype=float);d=x.sum(1).A1
    return x.multiply((1e6/d)[:,None]).tocsr()


def audit_results():
    m.check();from v62_conserved_response import data
    from run_v37_ensemble import checked_folder
    d=data();genes=d['genes'];targets=d['targets'];verified=0;manual=0
    for stage in ('diagnostic','development','confirmation'):
        for ctx in contexts():
            folder,_=checked_folder(ctx['family']);reference=sparse.load_npz(folder/'reference.npz');ref=normalized(reference);refmean=ref.mean(0).A1
            files=[p for p in (m.A/stage).glob(ctx['family']+'__*.json') if not p.name.endswith('__panel.json')]
            for path in files:
                r=read(path);raw=path.with_suffix('.counts.npz');arr=path.with_suffix('.npz')
                assert sha(raw)==r['counts_sha256'] and sha(arr)==r['arrays_sha256']
                c=sparse.load_npz(raw);assert c.shape==(r['cells'],len(genes)) and c.dtype.kind in 'iu' and c.data.min(initial=0)>=0
                actual=sparse.load_npz(folder/f"{r['case_index']}_actual.npz")
                if stage=='diagnostic':actual=actual[1::2]
                own=int(np.flatnonzero(genes==r['target'])[0]);axis=np.arange(len(genes))!=own
                pp=np.log1p(50000*c.sum(0).A1/c.sum());tt=np.log1p(50000*actual.sum(0).A1/actual.sum())
                np.testing.assert_allclose(float(np.sum((pp[axis]-tt[axis])**2)),r['profile_square_error'],rtol=1e-12,atol=1e-10)
                cn=normalized(c);an=normalized(actual);lfc=np.log2((cn.mean(0).A1+1e-9)/(refmean+1e-9));truth_lfc=np.log2((an.mean(0).A1+1e-9)/(refmean+1e-9));eligible=refmean>5
                with np.load(arr,allow_pickle=False) as z:
                    np.testing.assert_allclose(z['pred'],pp,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(z['lfc'],lfc,rtol=1e-10,atol=1e-10)
                    np.testing.assert_allclose(z['truth_lfc'],truth_lfc,rtol=1e-10,atol=1e-10)
                    p=z['p'].copy();q=z['q'].copy();cols=np.flatnonzero(eligible)
                    valid=eligible&axis;real=valid&(z['truth_q']<.05);adjudicable=valid&np.isfinite(truth_lfc)&(truth_lfc!=0)
                    called=adjudicable&(q<.05);correct=called&(np.sign(lfc)==np.sign(truth_lfc));budget=max(int(real.sum()),int(called.sum()))
                    union=real|(valid&(q<.05));expected_fid=float(correct.sum()/budget) if budget else None
                    expected_jac=float((real&(q<.05)).sum()/union.sum()) if union.any() else 1.
                    pool=np.flatnonzero(real&adjudicable);order=np.lexsort((genes[pool],-abs(lfc[pool]),p[pool],q[pool],q[pool]>=.05));pool=pool[order]
                    expected_reach=None
                    if len(pool):
                        precision=np.cumsum(np.sign(lfc[pool])==np.sign(truth_lfc[pool]))/np.arange(1,len(pool)+1);passing=np.flatnonzero(precision>=.9);expected_reach=float((passing[-1]+1)/int(real.sum())) if len(passing) else 0.
                    denominator=float(abs(truth_lfc[real]).sum());expected_nmae=float(abs(lfc[real]-truth_lfc[real]).sum()/denominator) if real.sum()>=10 and denominator>0 else None
                    for metric,value in [('fidelity',expected_fid),('reach',expected_reach),('jaccard',expected_jac),('nmae',expected_nmae)]:
                        if value is None:assert r[metric] is None,(path,metric,r[metric])
                        else:np.testing.assert_allclose(r[metric],value,rtol=1e-8,atol=1e-10,err_msg=str((path,metric)))
                    ordered=np.argsort(p[cols],kind='stable');rawq=p[cols][ordered]*len(cols)/np.arange(1,len(cols)+1);adjusted=np.minimum(1,np.minimum.accumulate(rawq[::-1])[::-1]);expected=np.empty(len(cols));expected[ordered]=adjusted
                    np.testing.assert_allclose(q[cols],expected,rtol=1e-12,atol=1e-12)
                    # A direct pooled-rank normal approximation, including ties.
                    choose=cols[np.linspace(0,len(cols)-1,8,dtype=int)]
                    for g in choose:
                        a=np.log1p(cn[:,g].toarray().ravel());b=np.log1p(ref[:,g].toarray().ravel());v=np.concatenate([a,b]);n=len(a);k=len(b);N=n+k
                        ranks=stats.rankdata(v);u=ranks[:n].sum()-n*(n+1)/2
                        _,sizes=np.unique(v,return_counts=True);variance=n*k/12*((N+1)-np.sum(sizes**3-sizes)/(N*(N-1)))
                        probability=1. if variance<=0 else min(1.,2*stats.norm.sf(abs(u-n*k/2)/np.sqrt(variance)))
                        np.testing.assert_allclose(probability,p[g],rtol=1e-9,atol=1e-12);manual+=1
                verified+=1
            del reference;gc.collect()
    # Trace checks confirm that no recipient study contributes a shape descriptor.
    for ctx in contexts():
        delta,trace=m.pooled((m.study(ctx['family']),))
        assert all(m.study(name)!=m.study(ctx['family']) for name in trace['training_contexts'])
        assert np.isfinite(delta).all() and np.all(delta[:,np.isin(genes,targets)]==0)
    write(m.A/'results_audit.json',dict(passed=True,count_populations=verified,manual_rank_coordinates=manual,
       checks=['count hashes/integer axes','independent normalized predicted/truth LFC and bulk profiles','independent four DE metrics','independent BH','manual tied ranks','whole-study shape traces','panel-gene abstention'],oracle_cases_not_deployable=True,D3_opened=False))
    print('V69 count result audit passed',verified,'populations',manual,'manual rank coordinates',flush=True)


def csr_group(handle,start,stop,ptr):
    lo,hi=int(ptr[start]),int(ptr[stop]);return sparse.csr_matrix((handle['X/data'][lo:hi],handle['X/indices'][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,18533))


def audit_payload():
    manifest,genes,targets,model,arrays=load();ready=read(m.A/'submission_ready.json');assert sha(m.OUT/'prediction.vcc')==ready['sha256']
    assert read(m.A/'results_audit.json')['passed'];provenance=pd.read_csv(m.A/'prediction_control_rows.csv');assert len(provenance)==360000
    scanned=0;replayed=0;compared=0;unchanged=0
    with tempfile.TemporaryDirectory(prefix='v69_payload_audit_') as tmp:
        payload=Path(tmp)/'prediction.h5ad'
        result=subprocess.run([str(m.ROOT/'.competition_env/Scripts/python.exe'),str(m.ROOT/'unpack_v24_official_payload.py'),str(m.OUT/'prediction.vcc'),str(payload)],capture_output=True,text=True,encoding='utf-8',errors='replace')
        assert result.returncode==0,result.stderr[-2000:];meta=json.loads(result.stdout);assert meta==ready['package_meta']
        obj=ad.read_h5ad(payload,backed='r');np.testing.assert_array_equal(obj.var_names,genes);np.testing.assert_array_equal(obj.obs.context.astype(str),np.repeat(manifest['contexts'],120000));np.testing.assert_array_equal(obj.obs.target_gene.astype(str),np.tile(np.repeat(targets,400),3));assert obj.obs_names.is_unique;obj.file.close()
        with h5py.File(payload,'r') as packed,h5py.File(m.OUT/'prediction.h5ad','r') as generated,h5py.File(m.ROOT/'experiments/v66_assay_strength/prediction.h5ad','r') as parent:
            ptr=packed['X/indptr'][:];pptr=parent['X/indptr'][:];np.testing.assert_array_equal(ptr,generated['X/indptr'][:])
            assert len(ptr)==360001 and ptr[-1]==meta['nnz']
            for start in range(0,int(ptr[-1]),4000000):
                stop=min(start+4000000,int(ptr[-1]))
                for k in ('data','indices'):np.testing.assert_array_equal(packed['X/'+k][start:stop],generated['X/'+k][start:stop])
                compared+=stop-start
            for ci,context in enumerate(manifest['contexts']):
                import build_submission as b
                x,guides,library,*_=b.load_control(m.ROOT,context,genes,manifest);source=ad.read_h5ad(m.ROOT/f'context_{context}.h5ad',backed='r');ids=source.obs_names.to_numpy();source.file.close()
                prov=provenance.iloc[ci*120000:(ci+1)*120000];rows=prov.source_control_row.to_numpy(int)
                np.testing.assert_array_equal(prov.source_control_cell_id.to_numpy(str),ids[rows]);np.testing.assert_array_equal(prov.ntc_id.to_numpy(str),guides[rows]);np.testing.assert_array_equal(prov.original_total_counts.to_numpy(int),library[rows])
                params=dict(p=arrays['bulk'][ci],alpha=arrays['alpha'][ci])
                for ti,target in enumerate(targets):
                    start=(ci*300+ti)*400;actual=csr_group(packed,start,start+400,ptr);local=prov.iloc[ti*400:(ti+1)*400];chosen=local.source_control_row.to_numpy(int);assert len(np.unique(chosen))==400
                    assert actual.data.min(initial=1)>0 and np.isfinite(actual.data).all() and (actual.data==np.floor(actual.data)).all();totals=actual.sum(1).A1;np.testing.assert_array_equal(totals,local.predicted_total_counts);assert totals.min()>0 and totals.max()<=1000000
                    seed=SEED+ci*100000+ti;old=csr_group(parent,start,start+400,pptr)
                    if model['candidate']=='hybrid25':
                        changed=np.random.default_rng(seed+690069).permutation(400)[:100];keep=np.setdiff1d(np.arange(400),changed);assert (actual[keep]-old[keep]).nnz==0;unchanged+=len(keep)
                    if ti in (0,1,149,299):
                        own=int(np.flatnonzero(genes==target)[0]);delta=None if model['candidate']=='mean_bridge' else arrays['delta'][ti]
                        expected,_=m.bridge(x[chosen],arrays['effects'][ci,ti],own,seed+690000,params,delta)
                        if model['candidate']=='hybrid25':
                            from run_v69_experiment import hybrid
                            expected=hybrid(old,expected,seed)
                        assert (actual-expected).nnz==0;replayed+=1
                    scanned+=400
                del x;gc.collect();print('V69 payload audit passed context',context,flush=True)
    assert scanned==360000 and replayed==12 and compared==ready['package_meta']['nnz']
    write(m.A/'independent_validation.json',dict(passed=True,verified_rows=scanned,verified_groups=900,packaged_nonzero_values_checked=compared,exact_replayed_groups=replayed,unchanged_parent_rows_verified=unchanged,package_sha256=ready['sha256'],prediction_sha256=sha(m.OUT/'prediction.h5ad'),model_sha256=sha(m.A/'reference_model.json'),validator_sha256=sha(Path(__file__)),all_counts_integer_and_valid=True,all_parent_control_identities_checked=True,variable_library_sizes_allowed=True,no_hidden_accuracy_claim=True))
    print('V69 complete payload validation passed',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['results','payload']);a=p.parse_args()
    with threadpool_limits(limits=2):(audit_results if a.stage=='results' else audit_payload)()
