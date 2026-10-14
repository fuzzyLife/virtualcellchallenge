"""Complete reserved prompts from frozen choices without mutating V20 numerical code.

The original preparer's tuple/list equality assertion prevents a JSON roundtrip
on repeat preparation. This helper uses the already frozen, audited choices and
exact original source pairing rule; no fit or selection parameters are changed.
"""
import json,gc
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
from v20_methods import ROOT,A,R,FAMILIES,align_counts,load_plan
from fit_v7_response import raw_inputs

def main():
 plan=load_plan();assert json.loads((A/'selection.json').read_text())['new_candidate_passes'];assert not json.loads((ROOT/'experiments/v21_stack_programs/analysis/selection.json').read_text())['new_candidate_passes']
 folder=A/'inputs';sources=folder/'sources';_,genes,targets=b.read_spec(ROOT);vocab=json.loads((R/'genes.json').read_text());contexts={c['family']:c for c in plan['contexts']};needed=set();files=[]
 for family in FAMILIES:
  rec=json.loads((folder/f'{family}.json').read_text())
  for c in contexts[family]['cells']:
   if c['split']!='confirmation':continue
   ti=int(np.flatnonzero(targets==c['target'])[0]);source=rec['source_choices'][ti]
   if source:
    assert source!=family;needed.add((source,ti));needed.add((source,rec['shuffled_target_indices'][ti]))
 for source in FAMILIES:
  requested=sorted(i for f,i in needed if f==source)
  if not requested:continue
  x,meta,measured,_=raw_inputs(source,genes);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);left=np.asarray(contexts[source]['left']);control={v:left[batches[left]==v] for v in sorted(set(batches[left]))}
  for ti in requested:
   target=str(targets[ti]);side=sources/f'{source}_{target}.json';rng=np.random.default_rng(np.random.SeedSequence([20020,FAMILIES.index(source),ti]));rows=rng.permutation(np.flatnonzero((labels==target)&np.isin(batches,list(control))))[:400];assert len(rows)>=40;pool={k:list(rng.permutation(v)) for k,v in control.items()};paired=[]
   for row in rows:
    key=batches[row]
    if not pool[key]:pool[key]=list(rng.permutation(control[key]))
    paired.append(pool[key].pop())
   paired=np.asarray(paired);assert np.array_equal(batches[rows],batches[paired]);records={}
   for kind,ids in [('treated',rows),('control',paired)]:
    p=sources/f'{source}_{target}_{kind}.npz';expected=align_counts(x[ids],genes,vocab)
    if p.exists():
     actual=sparse.load_npz(p);assert actual.shape==expected.shape and (actual-expected).nnz==0
    else:sparse.save_npz(p,expected)
    records[kind]={'rows':ids.tolist(),'unique_original_cells':len(set(ids)),'sha256':b.sha256(p),'full_library_mean':float(meta.total_counts.iloc[ids].mean())};files.append({'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(p)})
   record={'family':source,'target':target,'records':records,'vocabulary_measured_mask':np.isin(vocab,genes[measured]).tolist(),'batches_match_by_row':True,'source_input_rows':len(rows)}
   if side.exists():assert json.loads(side.read_text())==record
   else:b.write_json(side,record)
   files.append({'path':str(side.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(side)});print('Reserved prompt verified',source,target,len(rows),flush=True)
  del x;gc.collect()
 b.write_json(A/'confirmation_input_preparation.json',{'complete':True,'needed_source_target_pairs':len(needed),'source_target_pairs':[[f,int(i)] for f,i in sorted(needed)],'input_files':files,'plan_sha256':b.sha256(A/'plan.json'),'script_sha256':b.sha256(Path(__file__)),'query_treatments_used':False,'frozen_source_choices_reused':True,'original_raw_counts_and_same_batch_pairing_reconstructed':True,'reason':'Complete only the frozen reserved source prompts; original preparer repeats a tuple/list JSON equality assertion.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
