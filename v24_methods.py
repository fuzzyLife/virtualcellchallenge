"""V24: compatible CD4 condition means with donor-consistency confidence."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse,stats
from scipy.special import ndtr
import build_submission as b
from fit_v7_response import read_cache,raw_inputs
from v16_methods import source_arrays,emit as old_emit
from v8_methods import gated
from build_cd4_model import load_subset,eligible
from analyze_target_rank_correlations import pairwise_spearman
from v22_methods import write,digest,logcpm
ROOT=b.ROOT;OUT=ROOT/'experiments/v24_cd4_response';A=OUT/'analysis'
FAMILIES=['H1','HCT116','HEK293T'];NAMES=['v16','cdf25','cdf50','cdf50_de','shuffled25','shuffled50','shuffled50_de','amplitude'];ELIGIBLE=['cdf25','cdf50','cdf50_de']

def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V24 dependency changed: '+name
 for name,h in p['source_inputs'].items():assert b.sha256(ROOT/name)==h
 return p

def training_data(data,production=False):
 if production:return data
 keep=~((data['row_families']=='CD4')&(data['row_contexts']==3));result=dict(data)
 for k in ['truth','errors','valid','q','knockdown','row_families','row_contexts','row_targets']:result[k]=data[k][keep]
 return result

def raw(family,genes):
 if family in FAMILIES:return raw_inputs(family,genes)
 folder=ROOT/'experiments/v22_state_response/analysis/cd4_D1_selected' if family=='CD4_Rest_D1' else A/'holdout'
 rec=json.loads((folder/'extraction.json').read_text());x=sparse.load_npz(folder/'counts.npz');meta=pd.read_csv(folder/'cells.csv');axis=pd.read_csv(folder/'genes.csv');assert b.sha256(folder/'counts.npz')==rec['counts_sha256'] and b.sha256(folder/'cells.csv')==rec['cells_sha256'];np.testing.assert_array_equal(axis.gene_name.to_numpy(str),genes);meta['sample']=meta.lane_id.astype(str)
 return x,meta,axis.measured.to_numpy(bool),rec

def make_teacher(production=False):
 load_plan();name='teacher_all' if production else 'teacher_development';folder=A/'models';folder.mkdir(exist_ok=True);file=folder/(name+'.npz');side=file.with_suffix('.json')
 if side.exists():
  rec=json.loads(side.read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json')
  with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
 data=read_cache();genes=data['genes'];targets=data['targets'];cond=['Rest','Stim8hr','Stim48hr'] if production else ['Rest','Stim8hr'];directory=ROOT/'reference/cd4_probe';m,d=load_subset(directory,'merged');pairs=['CE0006864_CE0008162','CE0008678_CE0010866'];rep=[load_subset(directory,v) for v in pairs];lookups=[{(str(row.target_contrast),str(row.culture_condition)):i for i,row in f.iterrows()} for f,_ in rep];gene_lookup={str(g):i for i,g in enumerate(genes)};target_lookup={str(t):i for i,t in enumerate(targets)};destination=np.array([gene_lookup.get(str(g),-1) for g in d['genes']]);src=np.flatnonzero(destination>=0);dest=destination[src];rep_gene_maps=[{str(g):i for i,g in enumerate(v['genes'])} for _,v in rep];rep_axis=[np.array([mp.get(str(g),-1) for g in d['genes'][src]]) for mp in rep_gene_maps];shape=(len(cond),len(targets),len(genes));value=np.zeros(shape,np.float32);confidence=np.zeros(shape,np.float32);confidence_de=np.zeros(shape,np.float32);precision=np.zeros(shape,np.float32);sources=[];own=np.array([gene_lookup[t] for t in targets])
 for ri,row in m.iterrows():
  condition=str(row.culture_condition);target=str(row.target_contrast_gene_name)
  if condition not in cond or target not in target_lookup or not bool(eligible(m.iloc[[ri]])[0]):continue
  ci=cond.index(condition);ti=target_lookup[target];y=np.asarray(d['log_fc'][ri,src],float);se=np.asarray(d['lfcSE'][ri,src],float);q=np.asarray(d['adj_p_value'][ri,src],float);valid=np.isfinite(y)&np.isfinite(se)&(se>0)&np.isfinite(q);yy=np.clip(np.nan_to_num(y),-3,3);ss=np.maximum(np.nan_to_num(se,nan=1e6,posinf=1e6),.001);strong=(q<.05)&valid;direction=np.clip(2*ndtr(abs(yy)/ss)-1,0,1)*strong;equivalent=np.clip(1-(abs(yy)+1.96*ss)/.1,0,1)*valid;conf=.5*np.maximum(direction,equivalent);key=(str(row.target_contrast),condition);repvals=[]
  for k,(rm,rd) in enumerate(rep):
   arr=np.zeros(len(src));err=np.full(len(src),1e6);qq=np.ones(len(src));present=np.zeros(len(src),bool)
   if key in lookups[k]:
    rr=lookups[k][key];axis=rep_axis[k];ok=axis>=0;arr[ok]=rd['log_fc'][rr,axis[ok]];err[ok]=rd['lfcSE'][rr,axis[ok]];qq[ok]=rd['adj_p_value'][rr,axis[ok]];present=ok&np.isfinite(arr)&np.isfinite(err)&(err>0)&bool(eligible(rm.iloc[[rr]])[0])
   repvals.append((np.clip(np.nan_to_num(arr),-3,3),np.maximum(np.nan_to_num(err,nan=1e6),.001),qq,present))
  l,le,lq,lv=repvals[0];r,re,rq,rv=repvals[1];both=lv&rv;ls=(lq<.05)&lv;rs=(rq<.05)&rv;agree=(np.sign(l)==np.sign(r))&(np.sign(l)==np.sign(yy));corroborated=both&ls&rs&agree;donor_conf=np.minimum(np.clip(2*ndtr(abs(l)/le)-1,0,1),np.clip(2*ndtr(abs(r)/re)-1,0,1));conf[corroborated]=np.maximum(conf[corroborated],donor_conf[corroborated]);opposed=both&ls&rs&(np.sign(l)!=np.sign(r));conf[opposed]=0;null=np.minimum(np.clip(1-(abs(l)+1.96*le)/.1,0,1),np.clip(1-(abs(r)+1.96*re)/.1,0,1));conf[both]=np.maximum(conf[both],null[both]);conf[~valid]=0;conf_de=.5*direction;conf_de[corroborated]=np.maximum(conf_de[corroborated],donor_conf[corroborated]);conf_de[opposed]=0;conf_de[~valid]=0
  source_own=np.flatnonzero(d['genes']==target);kd=float(np.clip(1-np.exp2(min(float(d['log_fc'][ri,source_own[0]]),0)),0,.99)) if len(source_own) and np.isfinite(d['log_fc'][ri,source_own[0]]) else .9;gain=float(np.clip(.9/max(kd,.35),.8,1.8));yy=np.clip(yy*gain,-3,3);w=np.where(valid,1/np.maximum(ss**2,.0025),0);previous=precision[ci,ti,dest];total=previous+w;value[ci,ti,dest]=np.divide(value[ci,ti,dest]*previous+yy*w,total,out=np.zeros(len(dest)),where=total>0);confidence[ci,ti,dest]=np.divide(confidence[ci,ti,dest]*previous+conf*w,total,out=np.zeros(len(dest)),where=total>0);confidence_de[ci,ti,dest]=np.divide(confidence_de[ci,ti,dest]*previous+conf_de*w,total,out=np.zeros(len(dest)),where=total>0);precision[ci,ti,dest]=total;sources.append({'source_row':int(ri),'target':target,'condition':condition,'donor_corroborated_genes':int(corroborated.sum()),'credible_null_genes':int((equivalent>0).sum()),'opposed_donor_genes':int(opposed.sum())})
 for ci in range(len(cond)):value[ci,np.arange(300),own]=0;confidence[ci,np.arange(300),own]=0;confidence_de[ci,np.arange(300),own]=0
 ids=np.array([1,2,3][:len(cond)]);profiles=data['profiles'][ids];axis=np.flatnonzero(~np.isin(genes,targets)&np.isfinite(profiles).all(0)&(np.nanmax(profiles,axis=0)>5));rank=np.array([stats.rankdata(p[axis],method='average')/len(axis) for p in profiles]);state=axis[np.argsort(np.ptp(rank,axis=0),kind='stable')[-512:]];result={'genes':genes,'targets':targets,'value':value,'confidence':confidence,'confidence_de':confidence_de,'profiles':profiles,'conditions':np.array(cond),'state_axis':state,'own':own};np.savez_compressed(file,**result);write(side,{'sha256':b.sha256(file),'plan_sha256':b.sha256(A/'plan.json'),'training_conditions':cond,'withheld_response_condition':None if production else 'Stim48hr','donor_pair_ids':pairs,'supported_targets':int((confidence.sum(axis=(0,2))>0).sum()),'source_rows':sources,'scope':'Within-condition donor evidence; Rest development overlaps public summary biology and is calibration, not blind donor validation. Fresh-condition response rows excluded.'});print('V24 teacher fitted',name,len(sources),'source rows',flush=True);return result

def predict(data,teacher,query,excluded=(),production=False,baseline=None):
 bank=source_arrays(training_data(data,production),list(excluded));base=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32) if baseline is None else baseline.copy();axis=~np.isin(data['genes'],data['targets']);other=[]
 for ci,f in [(0,'K562'),(4,'H1'),(5,'HCT116'),(6,'HEK293T')]:
  if f not in excluded:other.append(float(pairwise_spearman(query[axis],data['profiles'][ci,axis],min_pairs=100)[0]))
 broad=np.array([pairwise_spearman(query[axis],p[axis],min_pairs=100)[0] for p in teacher['profiles']],float);gate=float(np.clip((np.nanmax(broad)-np.nanmax(other))/.1,0,1)) if np.isfinite(broad).any() and np.isfinite(other).any() else 0.;state=teacher['state_axis'];rho=np.array([pairwise_spearman(query[state],p[state],min_pairs=100)[0] for p in teacher['profiles']],float);rho=np.nan_to_num(rho,nan=-1);weight=np.exp((rho-rho.max())/.1);weight/=weight.sum();conf=np.einsum('c,ctg->tg',weight,teacher['confidence']);conf_de=np.einsum('c,ctg->tg',weight,teacher['confidence_de']);numerator=np.einsum('c,ctg,ctg->tg',weight,teacher['confidence'],teacher['value']);numerator_de=np.einsum('c,ctg,ctg->tg',weight,teacher['confidence_de'],teacher['value']);mu=np.divide(numerator,conf,out=np.zeros_like(numerator),where=conf>0);mu_de=np.divide(numerator_de,conf_de,out=np.zeros_like(numerator_de),where=conf_de>0);mu_de=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(mu_de,bank['own'])],np.float32);mu=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(mu,bank['own'])],np.float32);correction=gate*conf*np.clip(mu-base,-.5,.5);correction[np.arange(300),bank['own']]=0;available=np.flatnonzero(abs(correction).sum(1)>0);permutation=np.arange(300)
 if len(available):permutation[available]=np.random.default_rng(24024).permutation(available)
 shuffled=correction[permutation].copy();shuffled[np.arange(300),bank['own']]=0;norm=abs(correction).sum(1);othernorm=abs(shuffled).sum(1);result={'v16':base}
 for alpha in [.25,.5]:
  suffix=int(alpha*100);result['cdf'+str(suffix)]=np.clip(base+alpha*correction,-3,3).astype(np.float32);result['shuffled'+str(suffix)]=np.clip(base+alpha*shuffled,-3,3).astype(np.float32)
 correction_de=gate*conf_de*np.clip(mu_de-base,-.5,.5);correction_de[np.arange(300),bank['own']]=0;available_de=np.flatnonzero(abs(correction_de).sum(1)>0);perm_de=np.arange(300)
 if len(available_de):perm_de[available_de]=np.random.default_rng(24025).permutation(available_de)
 shuffled_de=correction_de[perm_de].copy();shuffled_de[np.arange(300),bank['own']]=0;result['cdf50_de']=np.clip(base+.5*correction_de,-3,3).astype(np.float32);result['shuffled50_de']=np.clip(base+.5*shuffled_de,-3,3).astype(np.float32)
 amplitude=np.divide(abs(result['cdf50']).sum(1),np.maximum(abs(base).sum(1),1e-12));result['amplitude']=(base*amplitude[:,None]).astype(np.float32)
 for v in result.values():v[np.arange(300),bank['own']]=0
 return result,{'lineage_compatibility_gate':gate,'CD4_broad_rho':broad.tolist(),'other_rho':other,'state_rho':rho.tolist(),'condition_weights':weight.tolist(),'changed_targets_at_alpha50':int((abs(result['cdf50']-base).sum(1)>0).sum())}

def emit(template,effect,own,seed,prepared,atlas):return old_emit(template,np.arange(template.shape[0]),effect,own,seed,'residual10',prepared,atlas)
