"""Reconstruct frozen Stack raw inputs and verify family exclusion/missing-gene masks."""
import json,gc
import numpy as np
from scipy import sparse
import build_submission as b
from v20_methods import ROOT,A,R,FAMILIES,load_plan,load_query
from fit_v7_response import raw_inputs,read_cache
from v19_methods import load_model as program_model,control_stats,residual_matrix
from v6_methods import moments
from threadpoolctl import threadpool_limits

def compare_projection(projected,original,genes,vocab):
 lookup={g:i for i,g in enumerate(genes)};cols=np.array([lookup.get(g,-1) for g in vocab]);present=cols>=0;assert projected[:,~present].nnz==0;assert (projected[:,present]!=original[:,cols[present]]).nnz==0

def main():
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];vocab=json.loads((R/'genes.json').read_text());assert len(vocab)==15012;assert b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256'];checks=[];source_sides=[ROOT/v['path'] for v in plan['input_files'] if '/inputs/sources/' in v['path'] and v['path'].endswith('.json')]
 for family in FAMILIES:
  x,meta,measured,_=raw_inputs(family,genes);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);ctx=next(c for c in plan['contexts'] if c['family']==family);query=load_query(family);ids=query['query_rows'];cached=query['cached_rows'];assert set(ids)==set(ctx['left']) and np.all(labels[ids]=='Non-Targeting') and set(cached)<=set(ids) and len(set(cached))==400
  compare_projection(sparse.load_npz(A/'inputs'/f'{family}_counts.npz'),x[cached],genes,vocab);np.testing.assert_array_equal(query['query_mean'],moments(x[ids])[1]*1e6);model=program_model(family);assert family not in model['training_families'];np.testing.assert_array_equal(query['positions'],model['positions']);np.testing.assert_array_equal(query['basis'],model['basis']);stat=control_stats(x,meta,ids);score=residual_matrix(stat,model['positions'])/model['scale']@model['basis'];np.testing.assert_allclose(score,query['query_scores'],atol=1e-12,rtol=1e-12)
  inputrec=json.loads((A/'inputs'/f'{family}.json').read_text());assert all(s!=family for s in inputrec['source_choices'] if s)
  for side in source_sides:
   rec=json.loads(side.read_text())
   if rec['family']!=family:continue
   trt=np.array(rec['records']['treated']['rows']);ntc=np.array(rec['records']['control']['rows']);assert len(set(trt))==len(trt)>=40;assert np.all(labels[trt]==rec['target']) and np.all(labels[ntc]=='Non-Targeting');assert np.array_equal(batches[trt],batches[ntc]) and set(ntc)<=set(ctx['left']);assert np.array_equal(rec['vocabulary_measured_mask'],np.isin(vocab,genes[measured]))
   for kind,rows in [('treated',trt),('control',ntc)]:compare_projection(sparse.load_npz(side.with_name(side.stem+'_'+kind+'.npz')),x[rows],genes,vocab)
  checks.append({'family':family,'query_controls':len(ids),'cached_distinct_query_controls':400,'query_treatment_values_used':False});del x,stat;gc.collect();print('Stack independent raw inputs/masks/control states verified',family,flush=True)
 b.write_json(A/'input_independent_audit.json',{'passed':True,'checks':checks,'source_prompt_pairs_reconstructed':len(source_sides),'all_counts_reconstructed_from_original_matrices':True,'source_batch_matching_verified':True,'missing_vocabulary_inputs_are_padding_not_measured_zero':True,'all_query_controls_and_program_states_verified':True,'pretraining_overlap_audited':False,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(ROOT/'audit_v20_inputs.py')})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
