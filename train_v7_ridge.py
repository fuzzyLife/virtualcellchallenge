"""Nested four-penalty Tahoe/rank ridge diagnostic and source-held residual baseline."""
import argparse,json,time,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import A,ROOT,read_cache
from v7_response_methods import FAMILIES,summary_score
from embedding_response_methods import kernels,predict_ridge
PENALTIES=[.01,.1,1.,10.];AMPLITUDES=[.35,.65,1.]

def rectangular_mask(data):
    mask=np.zeros_like(data['valid']);genes=data['genes'];panel=np.isin(genes,data['targets'])
    for ci in np.unique(data['row_contexts']):
        rows=np.flatnonzero(data['row_contexts']==ci)
        columns=data['valid'][rows].all(axis=0)&~panel
        mask[rows]=columns
    return mask

def estimate(data,train,targets,contexts,penalty):
    rows=np.flatnonzero(train);panel=np.array([int(np.flatnonzero(data['genes']==t)[0]) for t in data['targets']])
    embedding={'tahoe':(data['tahoe_vectors'][panel],data['tahoe_available'][panel])}
    k,cross=kernels(data['row_targets'][rows],data['row_contexts'][rows],targets,contexts,data['profiles'],embedding,'tahoe_rank')
    return predict_ridge(k,cross,data['truth'][rows],data['ridge_valid'][rows],data['row_families'][rows],penalty)

def main():
    data=read_cache();data['ridge_valid']=rectangular_mask(data);records=[];choices={};arrays={};splits=[]
    for outer in FAMILIES:
        inner='CD4' if outer!='CD4' else 'HEK293T';train=~np.isin(data['row_families'],[outer,inner]);validation=np.flatnonzero(data['row_families']==inner)
        local=[]
        for penalty in PENALTIES:
            pred,supported=estimate(data,train,data['row_targets'][validation],data['row_contexts'][validation],penalty)
            for amplitude in AMPLITUDES:
                scores=[]
                for i,ri in enumerate(validation):
                    valid=data['valid'][ri]&supported&~np.isin(data['genes'],data['targets'])
                    if valid.sum()<100:continue
                    scores.append(summary_score(pred[i]*amplitude,data['truth'][ri],data['errors'][ri],data['q'][ri],valid))
                f=pd.DataFrame(scores);row={'outer':str(outer),'inner':inner,'penalty':penalty,'amplitude':amplitude,'objective':float(f.selection_objective.mean()),'predicted_calls':float(f.predicted_calls.mean()),'fidelity_proxy':float(f.fidelity_proxy.mean())}
                records.append(row);local.append(row)
            print('Ridge inner',outer,penalty,flush=True)
        chosen=max(local,key=lambda x:x['objective']);choices[str(outer)]=chosen
        train=data['row_families']!=outer;test=np.flatnonzero(~train)
        for ci in np.unique(data['row_contexts'][test]):
            pred,supported=estimate(data,train,np.arange(300),np.full(300,ci),chosen['penalty'])
            arrays[f'{outer}_{ci}_ridge']=pred*chosen['amplitude'];arrays[f'{outer}_{ci}_supported']=supported
        splits.append({'outer':str(outer),'inner':inner,'inner_train_rows':np.flatnonzero(~np.isin(data['row_families'],[outer,inner])).tolist(),'validation_rows':validation.tolist(),'outer_train_rows':np.flatnonzero(train).tolist(),'outer_test_rows':test.tolist()})
        pd.DataFrame(records).to_csv(A/'ridge_inner_scores.csv',index=False)
        print('Ridge outer',outer,'selected',chosen['penalty'],chosen['amplitude'],flush=True)
    np.savez_compressed(A/'ridge_outer_predictions.npz',**arrays)
    penalty=float(pd.Series([c['penalty'] for c in choices.values()]).mode().iloc[0]);amplitude=float(np.median([c['amplitude'] for c in choices.values()]))
    _,genes,_=build.read_spec(ROOT);profiles=[];original=[]
    for context in ['A','B','C']:
        x,_,depth,_,_,_=build.load_control(ROOT,context,genes,__import__('json').loads((ROOT/'manifest.json').read_text()))
        profiles.append(np.asarray(x.multiply((1e6/depth)[:,None]).mean(axis=0)).ravel());del x;gc.collect()
        original.append({'context':context,'sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
    first=len(data['profiles']);data['profiles']=np.r_[data['profiles'],np.asarray(profiles)];train=np.ones(len(data['truth']),bool)
    pred,supported=estimate(data,train,np.tile(np.arange(300),3),np.repeat(np.arange(first,first+3),300),penalty)
    np.savez_compressed(A/'ridge_production_predictions.npz',predictions=(pred*amplitude).reshape(3,300,len(genes)),supported=supported)
    build.write_json(A/'ridge_manifest.json',{'outer_choices':choices,'production_penalty':penalty,'production_amplitude':amplitude,
        'source_cache_sha256':build.sha256(A/'source_cache.npz'),'outer_predictions_sha256':build.sha256(A/'ridge_outer_predictions.npz'),
        'production_predictions_sha256':build.sha256(A/'ridge_production_predictions.npz'),'script_sha256':build.sha256(ROOT/'train_v7_ridge.py'),
        'kernel_methods_sha256':build.sha256(ROOT/'embedding_response_methods.py'),'context_original_inputs':original,
        'mask':'Conservative intersection of observed genes across targets within each source context; all 300 panel output genes excluded only from this residual. Other responses remain supplied by pooling. No missing response labels imputed.',
        'selection':'Four penalties .01/.1/1/10 and amplitudes .35/.65/1; single inner whole-lineage selection excludes outer+inner; modal penalty and median amplitude for production.',
        'scope':'Public-source development regularization diagnostic; no hidden query treatment labels.'})
    build.write_json(A/'ridge_splits.json',splits)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-for-cache',action='store_true');a=p.parse_args()
    if a.wait_for_cache:
        while not (A/'source_provenance.json').exists():time.sleep(20)
    with threadpool_limits(limits=2):main()
