"""Prepare batch-paired public source prompts and controls-only Stack queries before freeze."""
from pathlib import Path
import json,gc
import numpy as np
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
import build_submission as b
from v20_methods import ROOT,A,R,FAMILIES,align_counts
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v19_methods import control_stats,residual_matrix,load_model as program_model

def main(stage='development'):
 A.mkdir(parents=True,exist_ok=True);folder=A/'inputs';folder.mkdir(exist_ok=True);sources=folder/'sources';sources.mkdir(exist_ok=True);data=read_cache();genes=data['genes'];targets=data['targets'];vocab=json.loads((R/'genes.json').read_text());parent=ROOT/'experiments/v19_covariance_forces/analysis';old=json.loads((parent/'plan.json').read_text());contexts={c['family']:c for c in old['contexts']};counts={};profiles={};masks={};metadata={}
 if stage=='confirmation':assert json.loads((A/'selection.json').read_text())['new_candidate_passes']
 for fi,family in enumerate(FAMILIES):
  x,meta,measured,rec=raw_inputs(family,genes);ids=np.array(contexts[family]['left']);profiles[family]=moments(x[ids])[1]*1e6;masks[family]=measured;metadata[family]=rec;valid_batch=set(meta['sample'].iloc[ids].astype(str));labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);counts[family]=np.array([np.count_nonzero((labels==t)&np.isin(batches,list(valid_batch))) for t in targets]);del x;gc.collect()
 needed=set()
 for fi,family in enumerate(FAMILIES):
  ctx=contexts[family];x,meta,measured,_=raw_inputs(family,genes);model=program_model(family);assert family not in model['training_families'];stat=control_stats(x,meta,ctx['left']);score=residual_matrix(stat,model['positions'])/model['scale']@model['basis'];cached=np.asarray(next(r['template'] for r in ctx['cells'] if r['split']=='development'));lookup={int(r):i for i,r in enumerate(ctx['left'])};cached_scores=score[[lookup[int(r)] for r in cached]];y=stat['y'][:,model['positions']].copy();y.data=np.expm1(y.data);y.data/=1+y.data;derivative=np.asarray(y.mean(0)).ravel()*np.log(2)
  qfile=folder/f'{family}_query.npz';matrix=folder/f'{family}_counts.npz'
  if stage=='development':
   np.savez_compressed(qfile,query_rows=np.asarray(ctx['left']),query_scores=score,cached_rows=cached,cached_scores=cached_scores,positions=model['positions'],basis=model['basis'],scale=model['scale'],derivative=derivative,query_mean=profiles[family],genes=genes,targets=targets,measured=measured);sparse.save_npz(matrix,align_counts(x[cached],genes,vocab))
  else:assert qfile.exists() and matrix.exists()
  similarities=[]
  for source in FAMILIES:
   if source==family:continue
   valid=measured&masks[source]&~np.isin(genes,targets)&((profiles[family]>0)|(profiles[source]>0));rho=float(stats.spearmanr(profiles[family][valid],profiles[source][valid]).statistic);similarities.append((source,rho,int(valid.sum())))
  similarities.sort(key=lambda v:(-v[1],v[0]));choices=[];shuffle=[];sample_sizes=[]
  for ti,target in enumerate(targets):
   available=[s for s,_,_ in similarities if counts[s][ti]>=40];source=available[0] if available else '';choices.append(source)
   if source:
    valid=np.flatnonzero(counts[source]>=40);bucket=int(np.digitize(counts[source][ti],[80,200]));group=valid[np.digitize(counts[source][valid],[80,200])==bucket]
    if len(group)>1:
     order=np.random.default_rng(20020+FAMILIES.index(source)*10+bucket).permutation(group);mapping=dict(zip(order,np.roll(order,1)));other=int(mapping[ti])
    else:
     choices_other=valid[valid!=ti];other=int(choices_other[np.argmin(abs(np.log(counts[source][choices_other]/counts[source][ti])))])
    shuffle.append(other);sample_sizes.append(int(min(400,counts[source][ti],counts[source][other])))
   else:shuffle.append(-1);sample_sizes.append(0)
  file=folder/f'{family}.json';record={'family':family,'query_arrays_sha256':b.sha256(qfile),'query_counts_sha256':b.sha256(matrix),'program_model_sha256':b.sha256(parent/'models'/f'{family}.npz'),'source_control_rank_correlations':similarities,'source_choices':choices,'shuffled_target_indices':shuffle,'matched_prompt_sizes':sample_sizes,'source_metadata':metadata,'query_treatments_used':False,'cached_query_unique_cells':len(set(cached))}
  if file.exists() and (A/'plan.json').exists():assert json.loads(file.read_text())==record
  else:b.write_json(file,record)
  for case in [r for r in ctx['cells'] if r['split']==stage]:
   ti=int(np.flatnonzero(targets==case['target'])[0]);source=choices[ti]
   if source:needed.add((source,ti));needed.add((source,shuffle[ti]))
  del x,stat;gc.collect();print('Stack query/input choices prepared',family,flush=True)
 for source in FAMILIES:
  requested=sorted(i for f,i in needed if f==source)
  if not requested:continue
  x,meta,measured,_=raw_inputs(source,genes);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);left=np.asarray(contexts[source]['left']);control={v:left[batches[left]==v] for v in sorted(set(batches[left]))}
  for ti in requested:
   target=str(targets[ti]);side=sources/f'{source}_{target}.json'
   if side.exists():continue
   rng=np.random.default_rng(np.random.SeedSequence([20020,FAMILIES.index(source),ti]));rows=np.flatnonzero((labels==target)&np.isin(batches,list(control)));rows=rng.permutation(rows)[:400];assert len(rows)>=40;pool={k:list(rng.permutation(v)) for k,v in control.items()};paired=[]
   for row in rows:
    key=batches[row]
    if not pool[key]:pool[key]=list(rng.permutation(control[key]))
    paired.append(pool[key].pop())
   paired=np.asarray(paired);assert np.array_equal(batches[rows],batches[paired]);records={}
   for kind,ids in [('treated',rows),('control',paired)]:
    p=sources/f'{source}_{target}_{kind}.npz';sparse.save_npz(p,align_counts(x[ids],genes,vocab));records[kind]={'rows':ids.tolist(),'unique_original_cells':len(set(ids)),'sha256':b.sha256(p),'full_library_mean':float(meta.total_counts.iloc[ids].mean())}
   b.write_json(side,{'family':source,'target':target,'records':records,'vocabulary_measured_mask':np.isin(vocab,genes[measured]).tolist(),'batches_match_by_row':True,'source_input_rows':len(rows)});print('Stack source paired',source,target,len(rows),flush=True)
  del x;gc.collect()
 b.write_json(A/f'{stage}_input_preparation.json',{'complete':True,'needed_source_target_pairs':len(needed),'source_target_pairs':[[f,int(i)] for f,i in sorted(needed)],'source_cache_sha256':b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'script_sha256':b.sha256(Path(__file__))})
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation'],nargs='?',default='development');args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage)
