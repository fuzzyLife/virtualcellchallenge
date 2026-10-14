"""Production inputs and bounded V20 mean responses, gated by reserved verification."""
from pathlib import Path
import json,gc
import numpy as np
from scipy import sparse,stats
import build_submission as b
from v20_methods import ROOT,A,R,FAMILIES,target_effects,load_plan
from fit_v7_response import raw_inputs
from v6_methods import moments
from v19_methods import load_model as program_model

PROD_DEPS=['v20_production.py','infer_v20_submission.py','fit_v20_submission.py','build_v20_submission.py','package_v20_submission.py','validate_v20_submission.py']

def approved():
 plan=load_plan()
 for name in ['selection.json','confirmation_result.json']:
  result=json.loads((A/name).read_text());assert result['new_candidate_passes'] and result['candidate']=='stack_mean';assert result['plan_sha256']==b.sha256(A/'plan.json')
 for name in ['confirmation_independent_audit.json','confirmation_source_input_independent_audit.json']:assert json.loads((A/name).read_text())['passed']
 return plan

def fast_align(x,genes,vocab):
 mapping={g:i for i,g in enumerate(vocab)};keep=np.array([i for i,g in enumerate(genes) if g in mapping]);z=x[:,keep].tocsr();indices=np.array([mapping[genes[i]] for i in keep]);out=sparse.csr_matrix((z.data.copy(),indices[z.indices],z.indptr.copy()),shape=(x.shape[0],len(vocab)));out.sort_indices();return out

def choices(query,measured,profiles,masks,counts,genes,targets,exclude=''):
 nonpanel=~np.isin(genes,targets);similarities=[]
 for source in FAMILIES:
  if source==exclude:continue
  valid=measured&masks[source]&nonpanel&((query>0)|(profiles[source]>0));similarities.append((source,float(stats.spearmanr(query[valid],profiles[source][valid]).statistic),int(valid.sum())))
 similarities.sort(key=lambda v:(-v[1],v[0]));selected=[];shuffled=[];sizes=[]
 for ti,target in enumerate(targets):
  available=[s for s,_,_ in similarities if counts[s][ti]>=40];source=available[0] if available else '';selected.append(source)
  if not source:shuffled.append(-1);sizes.append(0);continue
  valid=np.flatnonzero(counts[source]>=40);bucket=int(np.digitize(counts[source][ti],[80,200]));group=valid[np.digitize(counts[source][valid],[80,200])==bucket]
  if len(group)>1:
   order=np.random.default_rng(20020+FAMILIES.index(source)*10+bucket).permutation(group);other=int(dict(zip(order,np.roll(order,1)))[ti])
  else:
   rest=valid[valid!=ti];other=int(rest[np.argmin(abs(np.log(counts[source][rest]/counts[source][ti])))])
  shuffled.append(other);sizes.append(int(min(400,counts[source][ti],counts[source][other])))
 return {'source_control_rank_correlations':similarities,'source_choices':selected,'shuffled_target_indices':shuffled,'matched_prompt_sizes':sizes}

def prepare():
 plan=approved();folder=A/'production_inputs';folder.mkdir(exist_ok=True);sources=folder/'sources';sources.mkdir(exist_ok=True);assert not (A/'production_plan.json').exists();manifest,genes,targets=b.read_spec(ROOT);genes=np.asarray(genes,dtype=str);targets=np.asarray(targets,dtype=str);vocab=json.loads((R/'genes.json').read_text());contexts={c['family']:c for c in plan['contexts']};profiles={};masks={};counts={};source_records={};inputs=[];needed=set()
 for family in FAMILIES:
  x,meta,measured,record=raw_inputs(family,genes);left=np.asarray(contexts[family]['left']);profiles[family]=moments(x[left])[1]*1e6;masks[family]=measured;valid_batch=set(meta['sample'].iloc[left].astype(str));labels=meta.gene_target.to_numpy(str);batch=meta['sample'].to_numpy(str);counts[family]=np.array([np.count_nonzero((labels==t)&np.isin(batch,list(valid_batch))) for t in targets]);source_records[family]=record;del x;gc.collect()
 # Exact agreement with all 900 frozen public source/shuffle/size choices.
 for family in FAMILIES:
  expected=choices(profiles[family],masks[family],profiles,masks,counts,genes,targets,family);old=json.loads((A/'inputs'/(family+'.json')).read_text())
  for key,value in expected.items():assert json.loads(json.dumps(value))==old[key],key
 model=program_model('all')
 for ci,context in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6;rows=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([20261020,ci])));measured=np.ones(len(genes),bool);rec=choices(query,measured,profiles,masks,counts,genes,targets);qfile=folder/(context+'_query.npz');matrix=folder/(context+'_counts.npz');np.savez_compressed(qfile,query_mean=query,genes=genes,targets=targets,measured=measured,positions=model['positions'],basis=model['basis'],scale=model['scale'],cached_rows=rows);aligned=fast_align(x[rows],genes,vocab)
  from v20_methods import align_counts
  assert (aligned-align_counts(x[rows],genes,vocab)).nnz==0;sparse.save_npz(matrix,aligned);rec.update(context=context,query_arrays_sha256=b.sha256(qfile),query_counts_sha256=b.sha256(matrix),query_treatments_used=False,query_unique_cells=400,query_control_seed=[20261020,ci],original_input_sha256=b.sha256(ROOT/('context_'+context+'.h5ad')));b.write_json(folder/(context+'.json'),rec);inputs.append({'context':context,'file':'context_'+context+'.h5ad','sha256':rec['original_input_sha256']})
  for ti,source in enumerate(rec['source_choices']):
   if source:needed.add((source,ti))
  del x;gc.collect();print('Production query prepared',context,flush=True)
 for source in FAMILIES:
  requested=sorted(i for family,i in needed if family==source)
  if not requested:continue
  x,meta,measured,_=raw_inputs(source,genes);labels=meta.gene_target.to_numpy(str);batch=meta['sample'].to_numpy(str);left=np.asarray(contexts[source]['left']);control={v:left[batch[left]==v] for v in sorted(set(batch[left]))};observed=set(genes[measured]);vocab_mask=np.array([g in observed for g in vocab],bool)
  for offset,ti in enumerate(requested):
   target=str(targets[ti]);rng=np.random.default_rng(np.random.SeedSequence([20020,FAMILIES.index(source),ti]));rows=rng.permutation(np.flatnonzero((labels==target)&np.isin(batch,list(control))))[:400];assert len(rows)>=40;pool={k:list(rng.permutation(v)) for k,v in control.items()};paired=[]
   for row in rows:
    key=batch[row]
    if not pool[key]:pool[key]=list(rng.permutation(control[key]))
    paired.append(pool[key].pop())
   paired=np.asarray(paired);assert np.array_equal(batch[rows],batch[paired]);record={'family':source,'target':target,'source_input_rows':len(rows),'batches_match_by_row':True,'vocabulary_measured_mask':vocab_mask.tolist(),'records':{}}
   for kind,ids in [('treated',rows),('control',paired)]:
    path=sources/(source+'_'+target+'_'+kind+'.npz');sparse.save_npz(path,fast_align(x[ids],genes,vocab));record['records'][kind]={'rows':ids.tolist(),'unique_original_cells':len(set(ids)),'sha256':b.sha256(path)}
   b.write_json(sources/(source+'_'+target+'.json'),record)
   if (offset+1)%20==0:print('Production source prompts',source,offset+1,'/',len(requested),flush=True)
  del x;gc.collect()
 files=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(p)} for p in sorted(folder.rglob('*')) if p.is_file()]
 b.write_json(A/'production_plan.json',{'candidate':'stack_mean','development_plan_sha256':b.sha256(A/'plan.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'dependencies':{n:b.sha256(ROOT/n) for n in PROD_DEPS},'input_files':files,'context_original_inputs':inputs,'source_metadata':source_records,'all_frozen_public_source_choices_reproduced':True,'private_query_seed':20261020,'inference_seed_rule':'20261020 + context_index*10000 + target_index','scope':'Same validated bounded mean correction. All public sources available for ABC. No active variation or new tuning.'});print('Production inputs frozen',flush=True)

def verify_production():
 approved();plan=json.loads((A/'production_plan.json').read_text())
 for name,value in plan['dependencies'].items():assert b.sha256(ROOT/name)==value,'Production dependency changed: '+name
 for item in plan['input_files']:assert b.sha256(ROOT/item['path'])==item['sha256'],'Production input changed: '+item['path']
 return plan

if __name__=='__main__':
 from threadpoolctl import threadpool_limits
 with threadpool_limits(limits=2):prepare()
