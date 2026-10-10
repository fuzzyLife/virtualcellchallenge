"""Training-family-only response confidence and pathway-imputation screens."""
import json,gc
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from scipy import sparse
import build_submission as build
from fit_v7_response import read_cache
from v7_response_methods import FAMILIES,collapse,combine,adaptive_shrink,summary_score
from v11_methods import consensus_gate
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
NAMES=['v11','heterogeneity1','heterogeneity4','calibrated_gain','pathway_fill','shuffled_pathway','target_shuffle','train_centroid','zero']
def graph(targets,shuffle=False):
 members=pd.read_csv(ROOT/'analysis/post_v11_research/reactome_membership.csv');paths=pd.read_csv(ROOT/'analysis/post_v11_research/reactome_pathways.csv')
 small=set(paths.loc[(paths.all_genes>=3)&(paths.all_genes<=200),'pathway_id'])
 groups={t:set(members.loc[members.gene_name.eq(t)&members.pathway_id.isin(small),'pathway_id']) for t in targets}
 g=np.zeros((len(targets),len(targets)))
 for i,t in enumerate(targets):
  for j,u in enumerate(targets):
   if i!=j and groups[t]&groups[u]:g[i,j]=len(groups[t]&groups[u])/len(groups[t]|groups[u])
 # Relabel graph within annotation-degree bins, retaining edge structure and sparsity.
 if shuffle:
  degree=np.array([len(groups[t]) for t in targets]);perm=np.arange(len(targets));rng=np.random.default_rng(2026100933)
  bins=np.digitize(degree,[1,3,6,10,20,50])
  for b in np.unique(bins):ids=np.flatnonzero(bins==b);perm[ids]=rng.permutation(ids)
  g=g[perm][:,perm]
 return g

def pool(cube,hetero=.25,edges=None):
 y,var,valid,kd=[v.copy() for v in cube];gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 if edges is not None:
  for fi in range(len(y)):
   present=valid[fi].sum(axis=1)>100
   weights=edges*present[None,:];den=weights.sum(axis=1);n=(weights>0).sum(axis=1)
   w=np.divide(weights,den[:,None],out=np.zeros_like(weights),where=den[:,None]>0)
   ok=(~present)&(n>=2);assay=valid[fi].any(axis=0)
   for ti in np.flatnonzero(ok):
    cols=assay;mean=w[ti]@y[fi];unc=.5**2+(w[ti]@(var[fi]+(y[fi]-mean)**2))
    # Missing gene votes do not contribute to an inferred value or uncertainty.
    wg=w[ti,:,None]*valid[fi];d=wg.sum(axis=0);normalized=np.divide(wg,d,out=np.zeros_like(wg),where=d>0)
    inferred=(normalized*y[fi]).sum(axis=0);v=.25+(normalized*(var[fi]+(y[fi]-inferred)**2)).sum(axis=0)
    cols=cols&(d>0)&((valid[fi]&(weights[ti,:,None]>0)).sum(axis=0)>=2)
    y[fi,ti,cols]=inferred[cols];var[fi,ti,cols]=v[cols];valid[fi,ti,cols]=True
 mean,variance,ok=combine((y,var,valid,kd),np.array([2.,.25,2.,1.,1.]))
 if hetero!=.25:
  confidence=np.clip(.04/np.maximum(var,.0025),.1,4.);w=np.array([2.,.25,2.,1.,1.])[:,None,None]*confidence*valid;den=w.sum(axis=0)
  w=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0)
  h=(w*(y-mean)**2).sum(axis=0);variance=np.maximum(variance+(hetero-.25)*h,.0025)
 eb,_=adaptive_shrink(mean,variance,ok);gate,agreement,support=consensus_gate(y,var,valid,kd,[2.,.25,2.,1.,1.],mean,variance)
 return np.clip(eb+.2*gate*(mean-eb),-3,3),ok

def nested_gain(data,outer,reserved):
 ratios=[];counts=[]
 for inner in FAMILIES:
  if inner==outer:continue
  train=(data['row_families']!=outer)&(data['row_families']!=inner)&~np.isin(data['row_targets'],reserved)
  y,v,mask,kd=collapse(data,train);gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];v*=gain[:,:,None]**2
  mean,variance,ok=combine((y,v,mask,kd),np.array([2.,.25,2.,1.,1.]))
  # One coefficient, fit on other training families and non-confirmation targets only.
  rows=np.flatnonzero((data['row_families']==inner)&~np.isin(data['row_targets'],reserved));ti=data['row_targets'][rows]
  m=data['valid'][rows]&ok[ti];p=mean[ti];truth=data['truth'][rows]
  w=m/np.maximum(.25+data['errors'][rows]**2+variance[ti],.25)
  numerator=float((w*p*truth).sum());den=float((w*p*p).sum())
  if den>0:ratios.append(numerator/den);counts.append(int(m.sum()))
 return float(np.clip(np.mean(ratios),.25,1.25)),{'inner_ratios':ratios,'inner_pairs':counts}

def main():
 planpath=A/'response_screen_plan.json';plan={'candidates':NAMES,'whole_family_exclusion':True,'nested_gain_bounds':[.25,1.25],
  'pathway_policy':'Fill missing family/target only with >=2 annotated neighbors, Reactome paths 3..200 genes; output gene mask requires >=2 measured neighbor votes; uncertainty floor .5 LFC.',
  'negative_controls':['target label permutation','graph node permutation within annotation-degree bins','training centroid','zero'],
  'promotion':'Count stage remains frozen V11 biology. Response changes can advance only with consistent held-family improvements; no confirmation counts used here.',
  'script_sha256':build.sha256(Path(__file__)),'annotations_sha256':build.sha256(ROOT/'analysis/post_v11_research/reactome_membership.csv')}
 if planpath.exists():assert json.loads(planpath.read_text())==plan
 else:build.write_json(planpath,plan)
 data=read_cache();genes=data['genes'];targets=data['targets'];p=json.loads((A/'validation_plan.json').read_text());reserved={int(np.flatnonzero(targets==c['target'])[0]) for ctx in p['contexts'] for c in ctx['cells'] if c['split']=='confirmation'}
 edges=graph(targets);shuffled=graph(targets,True);records=[];details={}
 # All five source families are held out; evaluation targets exclude reserved count confirmation.
 for family in FAMILIES:
  cube=collapse(data,data['row_families']!=family);base,ok=pool(cube);bank={'v11':base}
  for name,h in [('heterogeneity1',1.),('heterogeneity4',4.)]:bank[name]=pool(cube,h)[0]
  gain,detail=nested_gain(data,family,list(reserved));bank['calibrated_gain']=base*gain;details[str(family)]={'gain':gain,**detail}
  bank['pathway_fill']=pool(cube,edges=edges)[0];bank['shuffled_pathway']=pool(cube,edges=shuffled)[0]
  perm=np.random.default_rng(2026100937).permutation(len(targets));bank['target_shuffle']=base[perm]
  tr=(data['row_families']!=family)&~np.isin(data['row_targets'],list(reserved));v=data['valid'][tr];den=v.sum(axis=0)
  centroid=np.divide((data['truth'][tr]*v).sum(axis=0),den,out=np.zeros(len(genes)),where=den>0)
  bank['train_centroid']=np.broadcast_to(centroid,base.shape);bank['zero']=np.zeros_like(base)
  rows=np.flatnonzero((data['row_families']==family)&~np.isin(data['row_targets'],list(reserved)))
  for ri in rows:
   ti=data['row_targets'][ri];mask=data['valid'][ri].copy();own=int(np.flatnonzero(genes==targets[ti])[0]);mask[own]=False
   for name,pred in bank.items():
    r=summary_score(pred[ti],data['truth'][ri],data['errors'][ri],data['q'][ri],mask)
    records.append({'family':str(family),'context':int(data['row_contexts'][ri]),'target':str(targets[ti]),'candidate':name,**r})
  pd.DataFrame(records).to_csv(A/'response_screen_scores.csv',index=False)
  print('Response screen completed',family,'nested gain',gain,flush=True);del cube,bank;gc.collect()
 df=pd.DataFrame(records);fields=['mse_ratio','fidelity_proxy','jaccard_proxy','sign_reference_DE','selection_objective']
 summary=df.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean();summary.to_csv(A/'response_screen_family.csv')
 summary.groupby('candidate').mean().to_csv(A/'response_screen_summary.csv')
 build.write_json(A/'response_screen_result.json',{'completed':True,'plan_sha256':build.sha256(planpath),'scores_sha256':build.sha256(A/'response_screen_scores.csv'),'calibration':details,
  'scope':'Existing source-summary normal-q proxies only; not competition count metrics, independent biological blind data, or hosted score estimates. Whole family excluded, reserved count-confirmation targets omitted from evaluation and gain fitting.'})
 print(summary.groupby('candidate').mean().to_string(),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
