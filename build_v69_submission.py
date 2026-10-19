"""Fit control parameters and stream a V69 artifact using frozen V66 means."""
from pathlib import Path
import argparse,gc,time
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,h5py,anndata as ad
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
import v69_count_bridge as m
from v57_pie import read,write,sha
from v6_methods import moments
SEED=20261011


def load():
    manifest,genes,targets=b.read_spec(m.ROOT);r=read(m.A/'reference_model.json')
    assert sha(m.A/'reference_model_arrays.npz')==r['arrays_sha256']
    assert sha(m.A/'plan.json')==r['plan_sha256']
    for name,h in r['dependencies'].items():assert sha(m.ROOT/name)==h,name
    with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
    np.testing.assert_array_equal(arrays['genes'],genes);np.testing.assert_array_equal(arrays['targets'],targets)
    return manifest,genes,targets,r,arrays


def fit(allow_pending_confirmation=False):
    m.check();selection=read(m.A/'selection.json');pending=not (m.A/'confirmation_assessment.json').exists()
    if pending:
        assert allow_pending_confirmation and selection['experimental_fallback'] and selection['selected']=='hybrid25'
        assessment=dict(selected=selection['selected'],release_supported=False)
    else:assessment=read(m.A/'confirmation_assessment.json')
    assert assessment['selected']==selection['selected'];assert not (m.A/'reference_model.json').exists()
    from build_v66_submission import load as parent_load
    manifest,genes,targets,parent,effects=parent_load();delta,trace=m.pooled(())
    alphas=[];bulk=[];cp=[];sources=[]
    for context in manifest['contexts']:
        x,guides,libs,*_=b.load_control(m.ROOT,context,genes,manifest);p=m.parameters(x)
        alphas.append(p['alpha']);bulk.append(p['p']);cp.append(moments(x)[1])
        sources.append(dict(context=context,file=f'context_{context}.h5ad',sha256=sha(m.ROOT/f'context_{context}.h5ad'),cells=x.shape[0]))
        del x;gc.collect()
    path=m.A/'reference_model_arrays.npz';np.savez_compressed(path,genes=np.asarray(genes,str),targets=np.asarray(targets,str),effects=effects,delta=delta,alpha=np.array(alphas),bulk=np.array(bulk),cp=np.array(cp))
    full_null=[read(m.A/'null'/f'{c}__bridge.json') for c in ('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2')]
    near_null=[read(m.A/'null'/f'{c}__near_zero_bridge.json') for c in ('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2')]
    null_ok=all(r['exact_count_identity'] and r['calls']<=5 for r in full_null) and all(r['calls']<=5 for r in near_null)
    supported=assessment['release_supported'] and null_ok and selection['selected']!='hybrid25'
    deps=['v69_count_bridge.py','build_v69_submission.py','run_v69_experiment.py']
    write(m.A/'reference_model.json',dict(candidate=selection['selected'],status='experimental_pending_confirmation' if pending else ('public_supported' if supported else 'experimental_failed_public_guards'),confirmation_pending=pending,
        available_public_guards_passed=supported,confirmed_improvement=False,arrays_sha256=sha(path),plan_sha256=sha(m.A/'plan.json'),
        selection_sha256=sha(m.A/'selection.json'),assessment_sha256=None if pending else sha(m.A/'confirmation_assessment.json'),source_trace=trace,
        parent_model_sha256=sha(m.ROOT/'experiments/v66_assay_strength/analysis/reference_model.json'),source_contexts=sources,
        dependencies={n:sha(m.ROOT/n) for n in deps},full_bridge_null_passed=null_ok,
        hybrid_null_identity=False,uploaded=False,D3_opened=False))
    print('V69 production controls/model fitted',selection['selected'],'supported',supported,flush=True)


def generate(preflight=False,workers=2):
    from v24_runtime import require_science
    require_science();manifest,genes,targets,record,model=load();choice=record['candidate'];records=[];provenance=[]
    parent=m.ROOT/'experiments/v66_assay_strength';parent_check=read(parent/'analysis/independent_validation.json')
    assert parent_check['passed'];parent_hash=sha(parent/'prediction.h5ad');assert parent_hash==parent_check['prediction_sha256']
    parent_prov=pd.read_csv(parent/'analysis/prediction_control_rows.csv')
    assert len(parent_prov)==360000
    writer=None;path=m.OUT/'prediction.h5ad'
    if not preflight:
        assert not path.exists() and not path.with_suffix('.h5ad.partial').exists() and not (m.OUT/'prediction.vcc').exists()
        obs=pd.DataFrame(dict(target_gene=pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),context=pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])),index=[f'v69_{i:06d}' for i in range(360000)])
        writer=b.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
    with h5py.File(parent/'prediction.h5ad','r') as src:
        ptr=src['X/indptr'][:]
        for ci,context in enumerate(manifest['contexts']):
            x,guides,library,*_=b.load_control(m.ROOT,context,genes,manifest)
            obj=ad.read_h5ad(m.ROOT/f'context_{context}.h5ad',backed='r');ids=obj.obs_names.to_numpy();obj.file.close()
            params=dict(p=model['bulk'][ci],alpha=model['alpha'][ci]);np.testing.assert_array_equal(params['p'],m.parameters(x)['p'])
            selected=['ABCD1','CHKB','VIM','FXYD5','HDAC8','NFE2L1','ACLY','ADNP']
            jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in selected]
            def compute(job):
                ti,target=job;start=(ci*300+ti)*400;prov=parent_prov.iloc[start:start+400];rows=prov.source_control_row.to_numpy(int)
                assert len(np.unique(rows))==400;np.testing.assert_array_equal(prov.source_control_cell_id.to_numpy(str),ids[rows]);np.testing.assert_array_equal(prov.original_total_counts.to_numpy(int),library[rows])
                seed=SEED+ci*100000+ti;own=int(np.flatnonzero(genes==target)[0]);template=x[rows]
                delta=None if choice=='mean_bridge' else model['delta'][ti]
                changed=400
                if choice=='hybrid25':
                    lo,hi=int(ptr[start]),int(ptr[start+400]);old=sparse.csr_matrix((src['X/data'][lo:hi],src['X/indices'][lo:hi],ptr[start:start+401]-lo),shape=(400,len(genes)))
                    selected_rows=np.random.default_rng(seed+690069).permutation(400)[:100]
                    generated,diag=m.bridge(template[selected_rows],model['effects'][ci,ti],own,seed+690000,params,delta)
                    new=old.tolil(copy=True);new[selected_rows]=generated;new=new.tocsr();changed=100
                else:new,diag=m.bridge(template,model['effects'][ci,ti],own,seed+690000,params,delta)
                totals=new.sum(1).A1;assert np.isfinite(new.data).all() and (new.data>0).all() and totals.min()>0 and totals.max()<=1000000
                mean=(new.multiply(1/totals[:,None])).mean(0).A1;bulk=new.sum(0).A1/totals.sum()
                factor=np.exp2(model['effects'][ci,ti]);factor[own]=.1;wanted=params['p']*factor;wanted/=wanted.sum()
                wanted_cp=model['cp'][ci]*factor;wanted_cp/=wanted_cp.sum()
                row=dict(mean_probability_L1=float(abs(mean-wanted_cp).sum()),context=context,target_gene=str(target),target_index=ti,nnz=int(new.nnz),bridge_cells=changed,mean_library=float(totals.mean()),maximum_library=int(totals.max()),mean_detected_genes=float(new.getnnz(1).mean()),bulk_probability_L1=float(abs(bulk-wanted).sum()),**{'sampled_bridge_'+k:v for k,v in diag.items()})
                pp=pd.DataFrame(dict(context=context,target_gene=str(target),prediction_cell_offset=np.arange(400),source_control_row=rows,source_control_cell_id=ids[rows],ntc_id=guides[rows],original_total_counts=library[rows],predicted_total_counts=totals))
                return new,row,pp
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
                for i in range(len(jobs)):
                    count,r,prov=pending.pop(i).result();records.append(r)
                    if writer is not None:writer.append(count);provenance.append(prov)
                    if i+workers<len(jobs):pending[i+workers]=pool.submit(compute,jobs[i+workers])
                    if preflight or (i+1)%10==0:print('V69 generation',context,i+1,'/',len(jobs),flush=True)
                    del count
            del x;gc.collect()
    pd.DataFrame(records).to_csv(m.A/('preflight_groups.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
    if preflight:
        write(m.A/'preflight.json',dict(passed=True,groups=len(records),candidate=choice));return
    writer.close();path.with_suffix('.h5ad.partial').replace(path)
    pd.concat(provenance,ignore_index=True).to_csv(m.A/'prediction_control_rows.csv',index=False)
    write(m.A/'run_parameters.json',dict(seed=SEED,candidate=choice,model_sha256=sha(m.A/'reference_model.json'),parent_prediction_sha256=parent_hash,generator_sha256=sha(Path(__file__)),source_row_provenance='Exact parent V66 control identities; quarter-row replacement if hybrid',groups=len(records),uploaded=False))
    write(m.A/'prediction_validation.json',b.validate_h5ad(path,genes,targets,manifest['contexts'],400))


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['fit','preflight','build']);p.add_argument('--workers',type=int,default=2);p.add_argument('--allow-pending-confirmation',action='store_true');a=p.parse_args()
    assert 1<=a.workers<=4
    with threadpool_limits(limits=2):
        if a.stage=='fit':fit(a.allow_pending_confirmation)
        else:generate(a.stage=='preflight',a.workers)
if __name__=='__main__':main()
