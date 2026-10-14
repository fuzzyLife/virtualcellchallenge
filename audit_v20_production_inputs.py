"""Independently compare every production Stack prompt with original raw cells."""
from pathlib import Path
import json,gc
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
from v20_production import ROOT,A,R,FAMILIES,verify_production
from fit_v7_response import raw_inputs

def compare_axis(original,stored,genes,vocab):
 mapping={str(g):i for i,g in enumerate(genes)};vi=np.array([i for i,g in enumerate(vocab) if g in mapping]);gi=np.array([mapping[vocab[i]] for i in vi]);assert stored.shape==(original.shape[0],len(vocab));assert (original[:,gi]-stored[:,vi]).nnz==0;missing=np.array([i for i,g in enumerate(vocab) if g not in mapping]);assert stored[:,missing].nnz==0;assert np.issubdtype(stored.dtype,np.integer) and (stored.data>=0).all()

def main():
 production=verify_production();folder=A/'production_inputs';manifest,genes,targets=b.read_spec(ROOT);vocab=json.loads((R/'genes.json').read_text());checks=[];plan=json.loads((A/'plan.json').read_text());contexts={c['family']:c for c in plan['contexts']};files=list((folder/'sources').glob('*.json'));scanned=0
 for context in manifest['contexts']:
  rec=json.loads((folder/(context+'.json')).read_text());x,guides,library,*_=b.load_control(ROOT,context,genes,manifest)
  with np.load(folder/(context+'_query.npz'),allow_pickle=False) as z:rows=z['cached_rows'].copy();query=z['query_mean'].copy();np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
  assert len(rows)==len(set(rows))==400;guide_counts=np.unique(guides[rows],return_counts=True)[1];assert len(guide_counts)==46 and guide_counts.min()==8 and guide_counts.max()==9;compare_axis(x[rows],sparse.load_npz(folder/(context+'_counts.npz')),genes,vocab);direct=np.asarray(x.astype(float).multiply((1e6/library)[:,None]).mean(0)).ravel();np.testing.assert_allclose(query,direct,rtol=1e-12,atol=1e-8);assert rec['query_unique_cells']==400 and not rec['query_treatments_used'];checks.append({'context':context,'query_cells_verified':400,'all_control_mean_CPM_independently_checked':True});del x;gc.collect();print('Production query raw-count audit',context,flush=True)
 for family in FAMILIES:
  selected=[p for p in files if p.stem.startswith(family+'_')]
  if not selected:continue
  x,meta,measured,_=raw_inputs(family,genes);labels=meta.gene_target.to_numpy(str);batch=meta['sample'].to_numpy(str);left=set(contexts[family]['left']);allowed_batches=set(batch[list(left)]);observed=set(genes[measured]);expected_mask=np.array([g in observed for g in vocab],bool)
  for offset,path in enumerate(selected):
   rec=json.loads(path.read_text());target=rec['target'];treated=np.asarray(rec['records']['treated']['rows']);control=np.asarray(rec['records']['control']['rows']);assert len(treated)>=40 and len(treated)==len(set(treated))==rec['source_input_rows'];assert np.all(labels[treated]==target) and set(control)<=left and np.all(labels[control]=='Non-Targeting');np.testing.assert_array_equal(batch[treated],batch[control]);np.testing.assert_array_equal(rec['vocabulary_measured_mask'],expected_mask);ti=int(np.flatnonzero(targets==target)[0]);expected=np.random.default_rng(np.random.SeedSequence([20020,FAMILIES.index(family),ti])).permutation(np.flatnonzero((labels==target)&np.isin(batch,list(allowed_batches))))[:400];np.testing.assert_array_equal(treated,expected)
   for kind,rows in [('treated',treated),('control',control)]:
    p=path.with_name(path.stem+'_'+kind+'.npz');assert b.sha256(p)==rec['records'][kind]['sha256'];compare_axis(x[rows],sparse.load_npz(p),genes,vocab);scanned+=len(rows)
   if (offset+1)%50==0:print('Production source raw-count audit',family,offset+1,'/',len(selected),flush=True)
  del x;gc.collect()
 b.write_json(A/'production_input_independent_audit.json',{'passed':True,'query_groups_verified':len(checks),'query_controls_verified':1200,'source_target_prompts_verified':len(files),'source_count_rows_verified_including_reused_controls':scanned,'all_vocab_count_values_and_missing_padding_compared_to_original_raw_counts':True,'all_source_target_labels_NTC_membership_unique_treatments_batch_pairing_and_seeds_verified':True,'production_plan_sha256':b.sha256(A/'production_plan.json'),'checks':checks,'auditor_sha256':b.sha256(Path(__file__))});print('All production inputs independently verified',len(files),'source prompts',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
