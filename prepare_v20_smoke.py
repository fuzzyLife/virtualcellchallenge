"""Build a bounded paired Stack smoke input from original public cells; no query treatments."""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
import build_submission as b
from fit_v7_response import raw_inputs,read_cache
ROOT=b.ROOT;A=ROOT/'experiments/v20_stack_transfer/analysis';R=ROOT/'reference/v20_stack'
def align(x,genes,vocab):
 mapping={g:i for i,g in enumerate(vocab)};keep=np.flatnonzero(np.isin(genes,vocab));z=x[:,keep].tocsr();indices=np.array([mapping[str(genes[i])] for i in keep]);out=sparse.csr_matrix((z.data.copy(),indices[z.indices],z.indptr.copy()),shape=(x.shape[0],len(vocab)));out.sort_indices();return out

def main():
 A.mkdir(parents=True,exist_ok=True);folder=A/'smoke_inputs';folder.mkdir(exist_ok=True);data=read_cache();genes=data['genes'];vocab=json.loads((R/'genes.json').read_text());plan=json.loads((ROOT/'experiments/v18_pie_expert/analysis/plan.json').read_text());ctx={c['family']:c for c in plan['contexts']};target='HDAC8';rng=np.random.default_rng(20020)
 x,meta,measured,_=raw_inputs('HEK293T',genes);labels=meta.gene_target.to_numpy(str);trt=np.flatnonzero(labels==target);assert len(trt)>=40;trt=rng.choice(trt,min(400,len(trt)),replace=False);ntc=np.asarray(ctx['HEK293T']['left']);ntc=rng.choice(ntc,400,replace=False);source_mask=np.isin(vocab,genes[measured]);records={}
 for name,rows in [('source_treated',trt),('source_control',ntc)]:
  z=align(x[rows],genes,vocab);p=folder/(name+'.npz');sparse.save_npz(p,z);records[name]={'family':'HEK293T','rows':rows.tolist(),'sha256':b.sha256(p),'cells':z.shape[0],'original_full_library_mean':float(meta.total_counts.iloc[rows].mean()),'vocabulary_library_mean':float(np.asarray(z.sum(1)).mean())}
 del x
 x,meta,measured,_=raw_inputs('H1',genes);case=next(r for r in ctx['H1']['cells'] if r['split']=='development' and r['target']==target);rows=np.asarray(case['template']);assert np.all(meta.gene_target.to_numpy(str)[rows]=='Non-Targeting');z=align(x[rows],genes,vocab);p=folder/'query_control.npz';sparse.save_npz(p,z);records['query_control']={'family':'H1','rows':rows.tolist(),'sha256':b.sha256(p),'cells':400,'original_full_library_mean':float(meta.total_counts.iloc[rows].mean()),'vocabulary_library_mean':float(np.asarray(z.sum(1)).mean())}
 supported=source_mask&np.isin(vocab,genes[measured]);b.write_json(folder/'inputs.json',{'target':target,'vocabulary':vocab,'supported_mask':supported.tolist(),'query_treatments_used':False,'records':records,'source_family_excluded_query':True,'code_sha256':b.sha256(Path(__file__))});print('Paired Stack raw-cell smoke inputs prepared',flush=True)
if __name__=='__main__':main()
