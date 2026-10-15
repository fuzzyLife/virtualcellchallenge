"""Attribute V26 null diagnostic errors to trained versus fallback coordinates."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v26_methods import ROOT,A,load_plan,load_decoder,geometry,raw,moments,state_values,emit,digest,logcpm,write
from v12_cpu_metrics import de
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b

def main():
 plan=load_plan();key='H1';ctx=next(c for c in plan['contexts'] if c['family']==key);case=ctx['null'];geo=geometry(key);model=load_decoder(key,'conditional');x,meta,measured,_=raw(key,geo['genes']);ids=case['fit_controls'];p,bulk=moments(x[ids])[1:3];state=state_values(x,meta,ids,model['state_axis'],model['state_basis'],model['state_scale']);lookup={v:i for i,v in enumerate(ids)};state=state[[lookup[v] for v in case['template']]];template=x[case['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);atlas,_=load_atlas_model();count,_=emit(template,np.zeros(len(geo['genes'])),None,case['seed'],'conditional',(libs,p,bulk),atlas,model,state);table=pd.read_csv(A/'null_H1_scores.csv');expected=table[table.candidate=='conditional'].iloc[0];assert digest(count)==expected.count_sha256;reference=x[case['reference']];score=de(count,reference,geo['genes']);calls=(score['q']<.05)&score['eligible'];assert calls.sum()==int(expected.null_DE_calls)
 training=sparse.load_npz(A/'training/H1_counts.npz')
 with np.load(A/'training/H1.npz',allow_pickle=False) as z:fit=z['fit']
 detected=np.zeros(len(geo['genes']),np.int64)
 for lo in range(0,len(fit),2048):
  rows=np.flatnonzero(fit[lo:lo+2048])+lo;detected+=np.asarray((training[rows]>0).sum(0)).ravel().astype(np.int64)
 panel=np.isin(geo['genes'],geo['targets']);category=np.where(panel,'panel_empirical_fallback',np.where(model['output_mask'],np.where(detected>=50,'trained_ge50_detections','trained_lt50_detections'),'other_empirical_fallback'));nc=logcpm(count);rc=logcpm(reference);nmean=np.asarray(nc.mean(0)).ravel();rmean=np.asarray(rc.mean(0)).ravel();nv=np.asarray(nc.power(2).mean(0)).ravel()-nmean*nmean;rv=np.asarray(rc.power(2).mean(0)).ravel()-rmean*rmean
 rows=pd.DataFrame({'gene_name':geo['genes'],'category':category,'fit_detected_cells':detected,'query_control_mean_CPM':p*1e6,'eligible':score['eligible'],'called':calls,'qvalue':score['q'],'LFC':score['lfc'],'generated_detection':np.asarray((count>0).mean(0)).ravel(),'reference_detection':np.asarray((reference>0).mean(0)).ravel(),'generated_log_variance':nv,'reference_log_variance':rv});rows[calls].sort_values('qvalue').to_csv(A/'H1_null_called_genes.csv',index=False);summary=rows.groupby('category').agg(genes=('gene_name','size'),eligible=('eligible','sum'),called=('called','sum')).reset_index();summary.to_csv(A/'H1_null_error_categories.csv',index=False);write(A/'H1_null_gene_diagnostic.json',{'plan_sha256':b.sha256(A/'plan.json'),'count_sha256':digest(count),'categories':summary.to_dict('records'),'called_genes_csv_sha256':b.sha256(A/'H1_null_called_genes.csv'),'diagnostic_sha256':b.sha256(Path(__file__)),'scope':'Post-rejection descriptive analysis of the already inspected H1 null. Does not change masks, sampling, weights or release thresholds; does not identify hidden competition errors.'});print(summary.to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
