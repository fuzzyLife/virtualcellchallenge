"""Nested, perturbation-conditioned neural response heads; no treated query inputs."""
from pathlib import Path
import argparse,json,time
import numpy as np,pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from neural_runtime import activate_local_runtime
activate_local_runtime()
import torch
from torch import nn
import build_submission as build
from fit_v7_response import A,read_cache,CONTEXTS
from v7_response_methods import FAMILIES,summary_score
SEED=20261008
ARCHITECTURES={'linear64':[64],'relu64':[64],'relu128_64':[128,64]}
def features(data,train,row_targets,row_contexts):
    genes=data['genes'];targets=data['targets'];gp={g:i for i,g in enumerate(genes)}
    panel=np.array([gp[t] for t in targets])
    unique=np.unique(data['row_targets'][train]);usable=unique[data['tahoe_available'][panel[unique]]]
    raw=data['tahoe_vectors'][panel].astype(float);center=raw[usable].mean(axis=0)
    _,_,vh=np.linalg.svd(raw[usable]-center,full_matrices=False)
    vectors=(raw-center)@vh[:16].T
    scale=np.maximum(vectors[usable].std(axis=0),.1);vectors/=scale
    vectors[~data['tahoe_available'][panel]]=0
    contexts=[]
    for profile in data['profiles']:
        values=profile[panel];valid=np.isfinite(values)
        rank=np.zeros(len(panel));rank[valid]=(rankdata(values[valid],method='average')-.5)/max(valid.sum(),1)
        contexts.append(np.r_[rank,valid.astype(float)])
    contexts=np.asarray(contexts)
    known=np.unique(data['row_contexts'][train]);cmean=contexts[known].mean(axis=0);cscale=np.maximum(contexts[known].std(axis=0),.1)
    contexts=(contexts-cmean)/cscale
    identity=np.eye(len(targets),dtype=np.float32)[row_targets]
    return np.c_[identity,vectors[row_targets],contexts[row_contexts]].astype(np.float32)
def make_model(size,architecture,outputs):
    dimensions=[size,*ARCHITECTURES[architecture],outputs];layers=[]
    for i,(left,right) in enumerate(zip(dimensions[:-1],dimensions[1:])):
        layers.append(nn.Linear(left,right))
        if architecture!='linear64' and i<len(dimensions)-2:layers.append(nn.ReLU())
    return nn.Sequential(*layers)
def loss_inputs(data,rows):
    y=data['truth'][rows].astype(np.float32);mask=data['valid'][rows].astype(np.float32)
    family=data['row_families'][rows];labels,counts=np.unique(family,return_counts=True)
    weight=np.array([len(rows)/len(labels)/counts[np.flatnonzero(labels==f)[0]] for f in family],dtype=np.float32)
    return y,mask,weight
def prediction_score(pred,data,rows):
    scores=[]
    for i,ri in enumerate(rows):
        valid=data['valid'][ri]
        if valid.sum()<100:continue
        scores.append(summary_score(pred[i],data['truth'][ri],data['errors'][ri],data['q'][ri],valid))
    frame=pd.DataFrame(scores)
    return float(frame.selection_objective.mean()),frame
def fit(data,train,architecture,epochs=80,validation=None):
    torch.manual_seed(SEED);torch.set_num_threads(2)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    rows=np.flatnonzero(train);target=data['row_targets'];context=data['row_contexts']
    x=features(data,train,target[rows],context[rows]);y,mask,weight=loss_inputs(data,rows)
    scale=np.sqrt(np.divide((y*y*mask).sum(axis=0),mask.sum(axis=0),out=np.zeros(y.shape[1]),where=mask.sum(axis=0)>0))
    scale=np.clip(scale,.25,2.).astype(np.float32)
    model=make_model(x.shape[1],architecture,y.shape[1]).to(device)
    xt=torch.tensor(x,device=device);yt=torch.tensor(y/scale,device=device)
    mt=torch.tensor(mask,device=device);wt=torch.tensor(weight,device=device)
    optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.01)
    rng=np.random.default_rng(SEED);records=[];best=(-np.inf,None)
    if validation is not None:
        vr=np.flatnonzero(validation)
        vx=torch.tensor(features(data,train,target[vr],context[vr]),device=device)
    for epoch in range(1,epochs+1):
        model.train()
        order=rng.permutation(len(rows))
        for start in range(0,len(rows),128):
            ix=order[start:start+128]
            prediction=model(xt[ix]);valid=mt[ix];weights=wt[ix,None]
            loss=(((prediction-yt[ix])**2)*valid*weights).sum()/torch.clamp((valid*weights).sum(),min=1)
            optimizer.zero_grad(set_to_none=True);loss.backward();optimizer.step()
        if validation is not None and epoch in [20,40,80]:
            model.eval()
            with torch.no_grad():pred=(model(vx).cpu().numpy()*scale)
            pred[:,~data['valid'][train].any(axis=0)]=0
            objective,frame=prediction_score(np.clip(pred,-3,3),data,vr)
            records.append({'architecture':architecture,'epoch':epoch,'objective':objective,'train_rows':len(rows),'validation_rows':len(vr)})
            if objective>best[0]:best=(objective,epoch)
    return model,scale,records,best
def predict(model,scale,data,train,targets,contexts):
    device=next(model.parameters()).device;f=features(data,train,targets,contexts);output=[]
    model.eval()
    with torch.no_grad():
        for start in range(0,len(f),128):output.append(model(torch.tensor(f[start:start+128],device=device)).cpu().numpy()*scale)
    result=np.clip(np.vstack(output),-3,3).astype(np.float32)
    result[:,~data['valid'][train].any(axis=0)]=0
    return result
def main():
    data=read_cache();records=[];scores=[];arrays={};choices={}
    for outer in FAMILIES:
        train=data['row_families']!=outer
        inner='HEK293T' if outer=='CD4' else 'CD4'
        fit_rows=train&(data['row_families']!=inner);validation=train&(data['row_families']==inner)
        assert not np.any(fit_rows&(data['row_families']==outer))
        candidates=[]
        for architecture in ARCHITECTURES:
            model,scale,r,best=fit(data,fit_rows,architecture,validation=validation)
            records.extend([{**record,'outer':str(outer),'inner':inner} for record in r])
            candidates.append((best[0],architecture,best[1]));del model
            print('Neural inner',outer,architecture,'objective',best[0],'epoch',best[1],flush=True)
        objective,architecture,epochs=max(candidates)
        choices[str(outer)]={'architecture':architecture,'epochs':epochs,'inner_family':inner,'inner_objective':objective}
        model,scale,_,_=fit(data,train,architecture,epochs)
        test=np.flatnonzero(~train)
        pred=predict(model,scale,data,train,data['row_targets'][test],data['row_contexts'][test])
        _,frame=prediction_score(pred,data,test);frame['family']=outer;frame['candidate']='nn_nested';scores.extend(frame.to_dict('records'))
        for ci in np.unique(data['row_contexts'][test]):
            arrays[f'{outer}_{ci}_nn_nested']=predict(model,scale,data,train,np.arange(len(data['targets'])),np.full(len(data['targets']),ci))
        pd.DataFrame(records).to_csv(A/'neural_inner_scores.csv',index=False);pd.DataFrame(scores).to_csv(A/'neural_outer_scores.csv',index=False)
        del model;torch.cuda.empty_cache()
        print('Neural outer complete',outer,architecture,flush=True)
    np.savez_compressed(A/'neural_outer_predictions.npz',**arrays)
    architecture=pd.Series([x['architecture'] for x in choices.values()]).mode().iloc[0]
    epochs=int(np.median([x['epochs'] for x in choices.values()]))
    train=np.ones(len(data['truth']),bool);model,scale,_,_=fit(data,train,architecture,epochs)
    weights={k:v.detach().cpu().numpy() for k,v in model.state_dict().items()}
    np.savez_compressed(A/'transition_weights.npz',scale=scale,**weights)
    build.write_json(A/'neural_manifest.json',{'outer_choices':choices,'production_architecture':architecture,'production_epochs':epochs,
        'source_cache_sha256':build.sha256(A/'source_cache.npz'),'prediction_sha256':build.sha256(A/'neural_outer_predictions.npz'),
        'weights_sha256':build.sha256(A/'transition_weights.npz'),'script_sha256':build.sha256(Path(__file__)),
        'torch_version':torch.__version__,'cuda':torch.cuda.is_available(),
        'inputs':'Target identity, frozen RNA-only target gene vectors projected on training targets, query-control panel ranks and measured-gene masks. No treated query.',
        'objective':'Masked source-family-balanced response loss; nested single-inner-lineage architecture/epoch selection by DE/error composite.',
        'scope':'Small supervised transition heads; not an Arc STATE foundation-model reproduction. Gene-feature transforms fitted on training contexts/targets.',
        'production_choice':'Modal inner-selected architecture and median epoch across outer folds; held-out lineage responses excluded from each inner fit/selection.'})
    print(pd.DataFrame(scores).groupby('family').mean(numeric_only=True).to_string(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-for-cache',action='store_true');a=p.parse_args()
    if a.wait_for_cache:
        while not (A/'source_provenance.json').exists():time.sleep(20)
    with threadpool_limits(limits=2):main()

