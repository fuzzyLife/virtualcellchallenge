"""V59: separate state abundance from within-state response; public data only."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(key,'2')
from pathlib import Path
from functools import lru_cache
import argparse,gc,tempfile
import numpy as np
import pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
from v57_pie import sha,read,write,ROOT,SOURCES
A=ROOT/'experiments/v59_state_composition/analysis'
CACHE=Path(tempfile.gettempdir())/'vcc_v59_state_composition'
GROUP={'H1':'H1','HCT116':'Orion','HEK293T':'Orion','CD4_Rest_D2':'CD4','CD4_Stim48_D2':'CD4'}
OUTERS=('H1','Orion','CD4');SEED=59059
ARMS=('v24','composition50','joint50','shuffled_joint50','pooled50')
CONFIG=dict(state_genes=512,latent=8,states=4,controls_per_context=512,state_pseudocount=.5,occupancy_prior_cells=32.,within_prior_cells=16.,log_tilt_bound=float(np.log(3)),within_lfc_bound=2.,correction_lfc_bound=1.,alpha=.5,count_seeds=[59059,159059],prediction_cells=400,maximum_reference_controls=2000)
def save(path,**arrays):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**arrays);return sha(path)
def load(path):
 with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
@lru_cache(None)
def plan():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,'V59 frozen dependency changed: '+n
 return p
def freeze():
 assert not (A/'plan.json').exists(),'Preserve prospective V59 plan'
 parent=read(ROOT/'experiments/v58_pie_holdstudy_guard/analysis/plan.json');inputs={}
 for s in SOURCES:
  f=ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(s+'.json');inputs[str(f.relative_to(ROOT))]=sha(f)
 deps=['v59_state_composition.py','test_v59_state_composition.py','evaluate_v59_counts.py','run_v59_experiment.py','v25_methods.py','fit_v7_response.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v57_pie.py']
 write(A/'plan.json',dict(version=59,seed=SEED,config=CONFIG,arms=ARMS,sources=SOURCES,groups=GROUP,inputs=inputs,dependencies={n:sha(ROOT/n) for n in deps},cohorts={s:dict(targets=parent['cohorts'][s]['targets'],count_pilot=parent['cohorts'][s]['count_pilot']) for s in SOURCES},
 state_model='512 training-only stable nonpanel genes; per-source/lane NTC centering/scaling; 8 PCA coordinates and 4 shared KMeans states. Complete outer study excluded from stable-gene selection, PCA, KMeans and source responses. Recipient NTCs only calibrate their own lane transform.',
 mechanism='Original V52 disjoint treatment/control views per lane. Smoothed state density ratios and within-state CPM log2FC; abundance shrink n/(n+32), within-state shrink state_n/(state_n+16). States from treated RNA can retain response. No observed individual counterfactual or birth/death claim.',
 oracle='Both view directions: fit on one treatment/control view, reconstruct the other using its controls, compare with other treatment. Composition, within, combined, direct mean and state-permuted composition. Diagnostic treatment evidence never enters deployable held-study inference.',
 transfer='Same-target source mechanisms; one total vote per study, reliability and nonpanel basal Spearman weighting. Measured masks prevent missing genes voting zero. Fixed primary joint50. Diagnostic composition50/pooled50. Within-source cyclic target derangement is negative control. All comparators use saved V58 baseline arrays; H1/Orion additionally exclude CD4 in these inherited comparator arrays, while held-CD4 comparator matches earlier tests.',
 release='Require fixed joint50 equal-study mean MSE gain >=1%, no study >2% worse, beat shuffled and pooled each study. Counts must improve expression MSE and two DE metrics, no CD4 FID/JAC regression >.005, beat shuffle. Public repeatedly examined biology; no hidden-score guarantee. D3 unopened. Failed candidates cannot be labeled improved or silently replace V24.',
 counts='44 fixed target-seed populations, common real rows/seeds across arms; unchanged V24 residual10 emitter for mean comparisons. New aggregate means permitted. State-cell diagnostic uses actual NTCs and conditional balanced count edits, without enforcing old means.',
 literature=[{'title':'Milo','url':'https://www.nature.com/articles/s41587-021-01033-z'},{'title':'muscat','url':'https://pubmed.ncbi.nlm.nih.gov/33257685/'}],D3_opened=False,automatically_upload=False))
 print('V59 frozen',sha(A/'plan.json'),flush=True)
@lru_cache(5)
def source_meta(source):
 f=ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(source+'.json');assert sha(f)==plan()['inputs'][str(f.relative_to(ROOT))];return read(f)
@lru_cache(1)
def axes():
 with np.load(ROOT/'experiments/v52_joint_mean_rank/analysis/labels/H1.npz',allow_pickle=False) as z:genes=z['genes'].copy();targets=z['targets'].copy()
 assert len(genes)==18533 and len(targets)==300
 return genes,targets
def prepare(source):
 from v25_methods import raw
 plan();folder=CACHE/'raw';folder.mkdir(parents=True,exist_ok=True);path=folder/(source+'.npz');receipt=A/'raw'/(source+'.json')
 if receipt.exists():assert sha(path)==read(receipt)['sha256'];return
 genes,targets=axes();x,meta,measured,r=raw(source,genes);labels=meta.gene_target.to_numpy(str);lane=meta['sample'].to_numpy(str);used=set(np.flatnonzero(labels=='Non-Targeting').tolist())
 for g in source_meta(source)['groups']:
  for q in g['lanes']:
   for role in ['treated_views','control_views']:
    for ids in q[role]:used.update(ids)
 ids=np.array(sorted(used),int);xx=x[ids].astype(float);depth=xx.sum(1).A1;assert (depth>0).all();xx=xx.multiply((1e6/depth)[:,None]).tocsr();sparse.save_npz(path,xx)
 side=folder/(source+'_meta.npz');save(side,rows=ids,label=labels[ids],lane=lane[ids],measured=measured,genes=genes,targets=targets,depth=depth)
 write(receipt,dict(sha256=sha(path),meta_sha256=sha(side),rows=len(ids),source_record=r,plan_sha256=sha(A/'plan.json'),normalization='Per-cell CPM over aligned gene axis as official DE; no rank pseudocount',D3_opened=False));print('Prepared',source,xx.shape,xx.nnz,flush=True);del x,xx;gc.collect()
def source_data(source):
 r=read(A/'raw'/(source+'.json'));path=CACHE/'raw'/(source+'.npz');assert sha(path)==r['sha256'];side=path.with_name(source+'_meta.npz');assert sha(side)==r['meta_sha256'];return sparse.load_npz(path),load(side)
def standardize(x,m,positions):
 result=np.zeros((x.shape[0],len(positions)),np.float32);transforms=[]
 for lane in sorted(set(m['lane'])):
  ids=np.flatnonzero(m['lane']==lane);ntc=ids[m['label'][ids]=='Non-Targeting'];assert len(ntc)>=2
  c=np.log1p(x[ntc][:,positions].toarray());mean=c.mean(0);sd=np.maximum(c.std(0),.5)
  for start in range(0,len(ids),1024):
   rows=ids[start:start+1024];result[rows]=np.clip((np.log1p(x[rows][:,positions].toarray())-mean)/sd,-5,5)
  transforms.append(dict(lane=str(lane),control_rows=m['rows'][ntc].tolist(),mean=mean.tolist(),sd=sd.tolist()))
 return result,transforms
def fit_states(outer):
 plan();path=CACHE/'states'/(outer+'.npz');rec=A/'states'/(outer+'.json')
 if rec.exists():assert sha(path)==read(rec)['sha256'];return
 genes,targets=axes();sources=[s for s in SOURCES if outer=='all' or GROUP[s]!=outer];stability=[];profiles=[];masks=[]
 for s in sources:
  f=ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(s+'.npz');assert sha(f)==source_meta(s)['sha256']
  with np.load(f,allow_pickle=False) as z:
   yi=z['mean'];valid=z['valid'];query=z['query'];energy=np.divide(np.sum(np.where(valid,abs(yi.mean(1)),0),axis=0),valid.sum(0),out=np.full(len(genes),np.inf),where=valid.sum(0)>0);stability.append(energy);profiles.append(query)
  _,m=source_data(s);masks.append(m['measured']);del m
 common=np.logical_and.reduce(masks)&~np.isin(genes,targets);common&=np.logical_and.reduce([np.isfinite(q)&(q>5) for q in profiles]);energy=np.mean(stability,axis=0);possible=np.flatnonzero(common&np.isfinite(energy));positions=possible[np.argsort(energy[possible],kind='stable')[:512]];assert len(positions)==512
 blocks=[];traces=[]
 for si,s in enumerate(sources):
  x,m=source_data(s);v,_=standardize(x,m,positions);ntc=np.flatnonzero(m['label']=='Non-Targeting');chosen=np.random.default_rng(SEED+si).choice(ntc,min(len(ntc),512),replace=False);blocks.append(v[chosen]);traces.append(dict(source=s,training_control_rows=m['rows'][chosen].tolist()));del x,v,m;gc.collect()
 matrix=np.concatenate(blocks);pca=PCA(n_components=8,svd_solver='full').fit(matrix);score=pca.transform(matrix);km=KMeans(n_clusters=4,n_init=10,random_state=SEED).fit(score)
 digest=save(path,positions=positions,pca_mean=pca.mean_,basis=pca.components_,centers=km.cluster_centers_,explained=pca.explained_variance_ratio_)
 write(rec,dict(sha256=digest,plan_sha256=sha(A/'plan.json'),outer=outer,training_sources=sources,source_controls=traces,stable_genes=genes[positions].tolist(),own_and_all_panel_genes_excluded=True,recipient_treatments_used=False,stable_gene_selection_uses_allowed_source_labels=True,explained_variance=pca.explained_variance_ratio_.tolist()));print('States fitted',outer,sources,flush=True)
def state_model(outer):
 p=CACHE/'states'/(outer+'.npz');assert sha(p)==read(A/'states'/(outer+'.json'))['sha256'];return load(p)
def assign(x,m,model):
 v,trans=standardize(x,m,model['positions']);z=(v-model['pca_mean'])@model['basis'].T;state=((z[:,None,:]-model['centers'][None,:,:])**2).sum(2).argmin(1);return state,trans
def group_summary(x,state,rows):
 rows=np.asarray(rows,int);n=np.bincount(state[rows],minlength=4);a=sparse.csr_matrix((np.ones(len(rows)),(state[rows],np.arange(len(rows)))),shape=(4,len(rows)));sums=(a@x[rows]).toarray();overall=x[rows].mean(0).A1;means=np.divide(sums,n[:,None],out=np.tile(overall,(4,1)),where=n[:,None]>0);return n.astype(float),means,overall
def estimate(nt,nc,tm,cm):
 tilt=np.log((nt+.5)/(nt.sum()+2))-np.log((nc+.5)/(nc.sum()+2));tilt-=tilt.mean();tilt=np.clip(tilt,-CONFIG['log_tilt_bound'],CONFIG['log_tilt_bound'])*nt.sum()/(nt.sum()+32)
 valid=(nt[:,None]>0)&(nc[:,None]>0)&(cm>5);within=np.clip(np.log2((tm+.1)/(cm+.1)),-2,2)*(nt/(nt+16))[:,None];return tilt,np.where(valid,within,0),valid
def reconstruct(nc,cm,tilt,within=None):
 weights=nc/max(nc.sum(),1)*np.exp(tilt);weights/=weights.sum();means=cm if within is None else cm*np.exp2(np.clip(within,-2,2));return weights@means
def logratio(a,b):return np.clip(np.log2((a+.1)/(b+.1)),-3,3)
def mechanisms(outer,source):
 plan();path=CACHE/'mechanisms'/outer/(source+'.npz');rec=A/'mechanisms'/outer/(source+'.json')
 if rec.exists():assert sha(path)==read(rec)['sha256'];return
 genes,targets=axes();x,m=source_data(source);state,trans=assign(x,m,state_model(outer));lookup=np.full(int(m['rows'].max())+1,-1,int);lookup[m['rows']]=np.arange(len(m['rows']))
 within=np.zeros((300,4,len(genes)),np.float32);valid=np.zeros_like(within,bool);tilt=np.zeros((300,4),np.float32);comp=np.zeros((300,len(genes)),np.float32);ncase=np.zeros(300,int);oracle=[];nonpanel=~np.isin(genes,targets)
 for case in source_meta(source)['groups']:
  ti=case['target_index'];target=case['target'];own=list(genes).index(target);numer=np.zeros((4,len(genes)));den=np.zeros_like(numer);lt=np.zeros(4);cc=np.zeros(len(genes));total=0
  for lane in case['lanes']:
   vv=[]
   for view in range(2):
    tr=lookup[np.array(lane['treated_views'][view],int)];cr=lookup[np.array(lane['control_views'][view],int)];assert (tr>=0).all() and (cr>=0).all() and not set(tr)&set(cr);assert np.all(m['label'][tr]==target) and np.all(m['label'][cr]=='Non-Targeting')
    nt,tm,tm0=group_summary(x,state,tr);nc,cm,cm0=group_summary(x,state,cr);l,w,ok=estimate(nt,nc,tm,cm);ok[:,~m['measured']]=False;ok[:,own]=False;w[~ok]=0;vv.append((nt,nc,tm,cm,tm0,cm0,l,w,ok));weight=len(tr);numer+=w*ok*weight;den+=ok*weight;lt+=l*weight;cc+=logratio(reconstruct(nc,cm,l),cm0)*weight;total+=weight
   for view in range(2):
    nt,nc,tm,cm,tm0,cm0,l,w,ok=vv[view];_,nc2,_,cm2,t2,c2,_,_,_=vv[1-view];axis=m['measured']&(c2>5)&nonpanel;axis[own]=False;truth=logratio(t2,c2)[axis]
    pred={'zero':np.zeros(len(genes)),'direct_mean':logratio(tm0,cm0),'composition':logratio(reconstruct(nc2,cm2,l),c2),'within':logratio(reconstruct(nc2,cm2,np.zeros(4),w),c2),'joint':logratio(reconstruct(nc2,cm2,l,w),c2),'state_permuted':logratio(reconstruct(nc2,cm2,np.roll(l,1)),c2)}
    for arm,y in pred.items():oracle.append(dict(source=source,target=target,lane=str(lane['lane']),fit_view=view,arm=arm,genes=int(axis.sum()),MSE=float(np.mean((y[axis]-truth)**2)),occupancy_view_difference=float(np.linalg.norm(vv[view][6]-vv[1-view][6]))))
  within[ti]=np.divide(numer,den,out=np.zeros_like(numer),where=den>0);valid[ti]=den>0;tilt[ti]=lt/max(total,1);comp[ti]=cc/max(total,1);ncase[ti]=total
 ntc=np.flatnonzero(m['label']=='Non-Targeting');nc,cm,ctrl=group_summary(x,state,ntc)
 digest=save(path,within=within,valid=valid,tilt=tilt,source_composition=comp,n=ncase,profile=ctrl,control_state_counts=nc,control_state_means=cm,state=state,rows=m['rows'],measured=m['measured'])
 op=A/'oracle'/outer/(source+'.csv');op.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(oracle).to_csv(op,index=False)
 write(rec,dict(sha256=digest,plan_sha256=sha(A/'plan.json'),source=source,outer=outer,targets=int((ncase>0).sum()),oracle_sha256=sha(op),state_transform_control_rows=trans,training_source=outer=='all' or GROUP[source]!=outer,diagnostic_only_when_recipient=True,D3_opened=False));print('Mechanisms',outer,source,int((ncase>0).sum()),flush=True);del x,within;gc.collect()
def mechanism(outer,source):
 p=CACHE/'mechanisms'/outer/(source+'.npz');assert sha(p)==read(A/'mechanisms'/outer/(source+'.json'))['sha256'];return load(p)
def spearman(a,b):
 good=np.isfinite(a)&np.isfinite(b)
 return float(np.corrcoef(rankdata(a[good],method='average'),rankdata(b[good],method='average'))[0,1]) if good.sum()>=100 else 0.
def transfer(recipient,outer=None):
 outer=GROUP[recipient] if outer is None else outer;plan();genes,targets=axes();recipient_all=mechanism(outer,recipient);q={k:recipient_all[k] for k in ['control_state_counts','control_state_means','profile','measured']};del recipient_all;nc=q['control_state_counts'];cm=q['control_state_means'];ctrl=q['profile'];source_names=[s for s in SOURCES if outer=='all' or GROUP[s]!=outer];donors={s:mechanism(outer,s) for s in source_names};axis=~np.isin(genes,targets);rhos={s:spearman(ctrl[axis],d['profile'][axis]) for s,d in donors.items()}
 parent=ROOT/'experiments/v58_pie_holdstudy_guard/analysis/arrays/inputs'/(recipient+'.npz')
 if not parent.exists():parent=Path(tempfile.gettempdir())/'vcc_v58_pie/inputs'/(recipient+'.npz')
 with np.load(parent,allow_pickle=False) as z:base=z['baseline'].copy()
 result={arm:base.copy() for arm in ARMS};tilts=np.zeros((300,4),np.float32);conditional=np.zeros((300,4,len(genes)),np.float32);available=np.zeros((300,len(genes)),bool);records=[]
 for ti,target in enumerate(targets):
  own=list(genes).index(target);values={}
  for shuffled in [False,True]:
   num=np.zeros((4,len(genes)));den=np.zeros_like(num);ln=np.zeros(4);ld=0.;sc=np.zeros(len(genes));sd=np.zeros(len(genes))
   for s,d in donors.items():
    j=ti
    if shuffled:
     ids=np.flatnonzero(d['n']>0)
     if ti not in ids:continue
     j=int(np.roll(ids,1)[np.flatnonzero(ids==ti)[0]])
    if d['n'][j]<=0:continue
    count=sum(GROUP[t]==GROUP[s] for t in source_names);weight=float(d['n'][j]/(d['n'][j]+64)*np.exp(2*rhos[s])/count);mask=d['valid'][j];num+=weight*d['within'][j]*mask;den+=weight*mask;ln+=weight*d['tilt'][j];ld+=weight;mm=d['measured'];sc+=weight*d['source_composition'][j]*mm;sd+=weight*mm
   ww=np.divide(num,den,out=np.zeros_like(num),where=den>0);ll=ln/max(ld,1e-12);sc=np.divide(sc,sd,out=np.zeros_like(sc),where=sd>0);mask=den.any(0)&q['measured'];mask[own]=False
   new=logratio(reconstruct(nc,cm,ll,ww),ctrl);pooled=logratio(reconstruct(nc,cm,np.zeros(4),ww),ctrl);comp=logratio(reconstruct(nc,cm,ll),ctrl)-sc;values[shuffled]=(new,pooled,comp,mask,ww,ll)
  new,pooled,comp,mask,ww,ll=values[False]
  for arm,y in [('joint50',new),('pooled50',pooled),('shuffled_joint50',values[True][0])]:
   support=values[True][3] if arm.startswith('shuffled') else mask;result[arm][ti,support]=np.clip(base[ti,support]+.5*np.clip(y[support]-base[ti,support],-1,1),-3,3)
  support=q['measured']&np.isfinite(comp);support[own]=False;result['composition50'][ti,support]=np.clip(base[ti,support]+.5*np.clip(comp[support],-1,1),-3,3)
  for v in result.values():v[ti,own]=0
  tilts[ti]=ll;conditional[ti]=ww;available[ti]=mask;records.append(dict(target=str(target),target_index=ti,supported_genes=int(mask.sum()),tilt=ll.tolist()))
 path=CACHE/'predictions'/(recipient+'.npz');digest=save(path,**result,tilt=tilts,within=conditional,available=available,genes=genes,targets=targets)
 write(A/'predictions'/(recipient+'.json'),dict(sha256=digest,outer=outer,training_sources=source_names,source_rank_correlations=rhos,records=records,plan_sha256=sha(A/'plan.json'),recipient_response_fields_used=False,D3_opened=False));print('Transferred',recipient,source_names,flush=True)
def evaluate():
 plan();rows=[];oracle=[]
 for s in SOURCES:
  z=load(ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(s+'.npz'));pr=load(CACHE/'predictions'/(s+'.npz'))
  for g in source_meta(s)['groups']:
   ti=g['target_index'];mask=z['valid'][ti]&~np.isin(z['genes'],z['targets']);truth=z['mean'][ti].mean(0)
   for arm in ARMS:
    err=pr[arm][ti,mask]-truth[mask];rows.append(dict(source=s,study=GROUP[s],target=g['target'],arm=arm,genes=int(mask.sum()),MSE=float(np.mean(err**2)),MAE=float(np.mean(abs(err))),supported_genes=int(pr['available'][ti,mask].sum())))
  oracle.append(pd.read_csv(A/'oracle'/GROUP[s]/(s+'.csv')))
 frame=pd.DataFrame(rows);frame.to_csv(A/'response_cases.csv',index=False);frame.groupby(['study','source','arm'])[['MSE','MAE']].mean().to_csv(A/'response_summary.csv');of=pd.concat(oracle,ignore_index=True);of.to_csv(A/'oracle_cases.csv',index=False);osum=of.groupby(['source','target','arm']).MSE.mean().groupby(['source','arm']).mean().unstack();osum.to_csv(A/'oracle_summary.csv')
 study=frame.groupby(['study','source','arm']).MSE.mean().groupby(['study','arm']).mean().unstack();ratio=study['joint50']/study['v24'];controls=(study['joint50']<study['shuffled_joint50'])&(study['joint50']<study['pooled50']);passed=bool(ratio.mean()<=.99 and (ratio<=1.02).all() and controls.all())
 write(A/'response_assessment.json',dict(passed=passed,joint50_MSE_ratio_by_study=ratio.to_dict(),equal_study_mean_ratio=float(ratio.mean()),beats_shuffled_and_pooled_each_study=controls.to_dict(),cases=len(frame)//len(ARMS),comparisons=len(frame),plan_sha256=sha(A/'plan.json'),response_sha256=sha(A/'response_cases.csv'),oracle_sha256=sha(A/'oracle_cases.csv'),count_tests_are_separate=True,D3_opened=False,private_score=None))
 print('Response MSE\n',study.to_string(),flush=True);print('Split oracle MSE\n',osum.to_string(),flush=True);print('Response gate:',passed,flush=True)
if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['freeze','prepare','states','mechanisms','transfer','evaluate']);q.add_argument('--source',choices=SOURCES);q.add_argument('--outer',choices=(*OUTERS,'all'));args=q.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='freeze':freeze()
  elif args.stage=='prepare':prepare(args.source)
  elif args.stage=='states':fit_states(args.outer)
  elif args.stage=='mechanisms':mechanisms(args.outer,args.source)
  elif args.stage=='transfer':transfer(args.source)
  else:evaluate()

