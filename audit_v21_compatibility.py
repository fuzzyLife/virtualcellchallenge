"""Independent dense query-state checks and parent-response invariance for V21."""
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
from v21_methods import ROOT,A,R,FAMILIES,load_plan,load_query,load_expert
from fit_v7_response import raw_inputs,read_cache
import v20_methods as parent

def main(stage):
 plan=load_plan();data=read_cache();vocab=json.loads((R/'genes.json').read_text());checks=[]
 for family in FAMILIES:
  q=load_query(family);file=A/'models'/f'{family}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json');assert family not in rec['training_families'];positions=q['positions'];assert np.isin(data['genes'][positions],vocab).all() and not np.isin(data['genes'][positions],data['targets']).any();np.testing.assert_allclose(q['basis'].T@q['basis'],np.eye(16),atol=1e-6)
  x,meta,measured,_=raw_inputs(family,data['genes']);ids=q['query_rows'];assert np.all(meta.gene_target.to_numpy(str)[ids]=='Non-Targeting');depth=meta.total_counts.to_numpy(float)[ids];batches=meta['sample'].to_numpy(str)[ids];names=sorted(set(batches));index=np.array([names.index(v) for v in batches]);y=np.log1p(x[ids][:,positions].toarray().astype(float)*1e4/depth[:,None]);means=np.array([y[index==i].mean(0) for i in range(len(names))]);z=np.log1p(depth);centers=np.array([z[index==i].mean() for i in range(len(names))]);zd=z-centers[index];beta=zd@y/max(float(zd@zd),len(ids)*1e-12);scores=(y-means[index]-zd[:,None]*beta)/q['scale']@q['basis'];np.testing.assert_allclose(scores,q['query_scores'],atol=1e-8,rtol=1e-8)
  new=load_expert(family,stage);oldfile=parent.A/'experts'/f'{family}_{stage}.npz';same=False
  if oldfile.exists():
   with np.load(oldfile,allow_pickle=False) as old:
    for key in ['lfc','shrink','raw','supported','shuffled_lfc','shuffled_shrink']:np.testing.assert_array_equal(new[key],old[key])
   same=True
  rows=np.flatnonzero(new['valid']);norm=np.linalg.norm(new['program_delta'][rows].reshape(len(rows),-1),axis=1);other=np.linalg.norm(new['shuffled_program_delta'][rows].reshape(len(rows),-1),axis=1);np.testing.assert_allclose(norm,other,rtol=1e-6);assert (norm>0).all();checks.append({'family':family,'all_program_genes_in_vocabulary':True,'query_scores_independently_recomputed':True,'mean_arrays_bitwise_match_parent':same,'active_targets':len(rows),'minimum_program_norm':float(norm.min()),'matched_shuffled_norm':True});del x;gc.collect();print('V21 compatibility/mean/dense state audit passed',family,flush=True)
 b.write_json(A/f'{stage}_compatibility_independent_audit.json',{'passed':True,'checks':checks,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(ROOT/'audit_v21_compatibility.py'),'pretraining_overlap_audited':False})
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage)
