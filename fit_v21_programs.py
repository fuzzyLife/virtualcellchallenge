"""Fit control-only shared programs on genes present in the pretrained Stack vocabulary."""
import json,gc
import numpy as np
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits
import build_submission as b
from v21_methods import ROOT,A,R,FAMILIES,load_plan
import v20_methods as parent
from v19_methods import control_stats,residual_matrix
from fit_v7_response import raw_inputs,read_cache

def main():
 plan=load_plan();data=read_cache();vocab=json.loads((R/'genes.json').read_text());contexts={c['family']:c for c in plan['contexts']};folder=A/'models';folder.mkdir(exist_ok=True);inputs=A/'inputs';inputs.mkdir(exist_ok=True)
 for held in FAMILIES+['all']:
  file=folder/f'{held}.npz'
  if file.exists():raise FileExistsError('Preserve fitted V21 programs')
  training=[f for f in FAMILIES if f!=held];stats={};masks=[]
  for family in training:
   x,meta,measured,_=raw_inputs(family,data['genes']);stats[family]=control_stats(x,meta,contexts[family]['left']);masks.append(measured);del x;gc.collect()
  variance=np.mean([stats[f]['variance'] for f in training],axis=0);mask=np.logical_and.reduce(masks)&np.isin(data['genes'],vocab)&~np.isin(data['genes'],data['targets'])&(variance>1e-4);possible=np.flatnonzero(mask);positions=possible[np.argsort(variance[possible],kind='stable')[-2048:]];assert len(positions)==2048;scale=np.sqrt(variance[positions]);matrix=np.concatenate([residual_matrix(stats[f],positions)/scale/np.sqrt(len(stats[f]['control_rows'])) for f in training]);_,values,components=randomized_svd(matrix,n_components=16,n_iter=3,random_state=19019);basis=components.T;del matrix;np.savez_compressed(file,positions=positions,genes=data['genes'][positions],scale=scale,basis=basis,singular_values=values,training_families=np.asarray(training,str));b.write_json(file.with_suffix('.json'),{'sha256':b.sha256(file),'plan_sha256':b.sha256(A/'plan.json'),'training_families':training,'all_program_genes_in_Stack':True,'all_panel_genes_excluded':True});del stats;gc.collect()
  if held!='all':
   old=parent.load_query(held);x,meta,_,_=raw_inputs(held,data['genes']);stat=control_stats(x,meta,old['query_rows']);score=residual_matrix(stat,positions)/scale@basis;lookup={int(r):i for i,r in enumerate(old['query_rows'])};cached_scores=score[[lookup[int(r)] for r in old['cached_rows']]];y=stat['y'][:,positions].copy();y.data=np.expm1(y.data);y.data/=1+y.data;derivative=np.asarray(y.mean(0)).ravel()*np.log(2);old.update(positions=positions,scale=scale,basis=basis,query_scores=score,cached_scores=cached_scores,derivative=derivative);qfile=inputs/f'{held}_query.npz';np.savez_compressed(qfile,**old);rec=json.loads((parent.A/'inputs'/f'{held}.json').read_text());rec.update(parent_query_arrays_sha256=rec['query_arrays_sha256'],query_arrays_sha256=b.sha256(qfile),plan_sha256=b.sha256(A/'plan.json'),program_model_sha256=b.sha256(file));b.write_json(inputs/f'{held}.json',rec);assert old['measured'][positions].mean()>=.8;del x,stat;gc.collect()
  print('V21 Stack-covered programs fitted',held,flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
