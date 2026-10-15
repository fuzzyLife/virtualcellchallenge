"""Descriptive localization of frozen V27 sampler artifacts, not new tuning."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v27_methods import ROOT,A,load_plan,write
from v26_methods import raw,moments
from fit_v7_response import read_cache
from v12_cpu_metrics import de
import build_submission as b

def main():
 plan=load_plan();data=read_cache();genes=data['genes'];case=plan['contexts'][0]['null_cases'][0];x,m,measured,_=raw('H1',genes);ref=x[case['reference']];base=sparse.load_npz(A/'sampling/raw_template.npz');labels=np.where(np.isin(genes,data['targets']),'panel','other');axis=measured&(ref.sum(0).A1>0);rows=[];tops=[];previous=None
 for name in pd.read_csv(A/'sampling_scores.csv').stage:
  c=sparse.load_npz(A/'sampling'/f'{name}.npz');v=c.toarray().astype(np.int64);r=ref.toarray().astype(np.int64);d=de(c,ref,genes);calls=(d['q']<.05)&d['eligible'];libs=v.sum(1);detect=(v>0).mean(0);real=(r>0).mean(0);total=v.sum(0);maxshare=np.divide(v.max(0),total,out=np.zeros(v.shape[1]),where=total>0)
  record={'stage':name,'signed_detection_error':float((detect-real)[d['eligible']].mean()),'mean_max_cell_share':float(maxshare[d['eligible']].mean()),'genes_ge_half_counts_one_cell':int(((maxshare>=.5)&d['eligible']).sum()),'changed_molecule_fraction_vs_raw':float(abs(v-base.toarray()).sum()/2/libs.sum()),'panel_DE_calls':int((calls&(labels=='panel')).sum()),'other_DE_calls':int((calls&(labels=='other')).sum())}
  if previous is not None:record.update({'newly_nonzero_entries_vs_previous':int(((v>0)&(previous==0)).sum()),'lost_nonzero_entries_vs_previous':int(((v==0)&(previous>0)).sum()),'changed_molecule_fraction_vs_previous':float(abs(v-previous).sum()/2/libs.sum())})
  rows.append(record)
  ids=np.flatnonzero(calls);ids=ids[np.argsort(d['q'][ids],kind='stable')[:50]]
  for j in ids:tops.append({'stage':name,'gene_name':genes[j],'category':labels[j],'q':float(d['q'][j]),'log2FC':float(d['lfc'][j]),'generated_detection':float(detect[j]),'reference_detection':float(real[j]),'largest_cell_count_share':float(maxshare[j])})
  previous=v
 pd.DataFrame(rows).to_csv(A/'sampling_artifact_diagnostic.csv',index=False);pd.DataFrame(tops).to_csv(A/'sampling_top_called_genes.csv',index=False);write(A/'sampling_artifact_diagnostic.json',{'plan_sha256':b.sha256(A/'plan.json'),'rows_sha256':b.sha256(A/'sampling_artifact_diagnostic.csv'),'genes_sha256':b.sha256(A/'sampling_top_called_genes.csv'),'diagnostician_sha256':b.sha256(Path(__file__)),'scope':'Spent H1 null, descriptive count changes and extreme cells. No new weights, parameters, selection or release thresholds.'});print(pd.DataFrame(rows).to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
