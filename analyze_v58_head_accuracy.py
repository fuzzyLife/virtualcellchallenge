"""Assess raw PIE DE/LFC heads on frozen real populations before generating cells."""
from functools import lru_cache
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from threadpoolctl import threadpool_limits
from v58_pie import *
from fit_v58_pie_decoder import inputs,predictions,predict
from v25_methods import raw
from v12_cpu_metrics import de
import build_v57_pie_evidence as evidence
from run_v58_cached_stage import ReadOnlyNumpy

def main():
 assert read(A/'response_evaluation.json')['complete']
 evidence.np=ReadOnlyNumpy();evidence.read=lru_cache(maxsize=None)(evidence.read)
 before=not (A/'counts').exists();rows=[];folder=CACHE/'head_truth';folder.mkdir(exist_ok=True)
 for source in SOURCES:
  z=inputs(source);pr=predictions(source);genes=z['genes'];x,meta,_,_=raw(source,genes);labels=meta.gene_target.to_numpy(str);ntc=np.flatnonzero(labels=='Non-Targeting')
  for target in plan()['cohorts'][source]['count_pilot']:
   ti=list(z['targets']).index(target);own=list(genes).index(target);mu,_,_=predict(source,ti);pos=pr['positions'];col=np.arange(len(pos));posmask=pos!=own
   for seed in (plan()['config']['count_seeds'] if source.startswith('CD4') else plan()['config']['count_seeds'][:1]):
    rng=np.random.default_rng(seed+ti);order=rng.permutation(ntc);template=order[:400];ref=order[400:2400];actual=rng.permutation(np.flatnonzero(labels==target))[:400]
    truth=de(x[actual],x[ref],genes);eligible=posmask&truth['eligible'][pos];called=eligible&(truth['q'][pos]<.05);n=int(called.sum());prob=pr['p_de'][ti,col][eligible];binary=(truth['q'][pos][eligible]<.05)
    r={'source':source,'target':target,'seed':seed,'seen':target in plan()['seen_targets'],'transfer':source.startswith('CD4'),'actual_cells':len(actual),'eligible_genes':int(eligible.sum()),'true_DE_genes':n,'DE_prevalence':float(binary.mean()),'raw_pDE_AP_eligible_only':float(average_precision_score(binary,prob)) if n else None}
    for name,values in [('raw_PIE_LFC',pr['lfc'][ti,col]),('V24_intended',z['baseline'][ti,pos]),('calibrated_intended',mu[pos])]:
     r[name+'_DE_MAE']=float(abs(values[called]-truth['lfc'][pos][called]).mean()) if n else None
     r[name+'_DE_direction']=float((np.sign(values[called])==np.sign(truth['lfc'][pos][called])).mean()) if n else None
    rows.append(r);file=folder/f'{source}_{target}_{seed}.npz';np.savez_compressed(file,q=truth['q'],lfc=truth['lfc'],eligible=truth['eligible'])
    write(A/'head_truth'/(file.stem+'.json'),{'sha256':sha(file),'plan_sha256':sha(A/'plan.json'),'actual_rows':actual.tolist(),'reference_rows':ref.tolist(),'template_rows':template.tolist()})
   print('Raw-head accuracy evaluated',source,target,flush=True)
 frame=pd.DataFrame(rows);frame.to_csv(A/'raw_head_cases.csv',index=False);frame.groupby(['source','seen']).mean(numeric_only=True).to_csv(A/'raw_head_summary.csv')
 write(A/'raw_head_evaluation.json',{'complete':True,'initial_run_before_count_generation':before,'cases':len(frame),'cases_sha256':sha(A/'raw_head_cases.csv'),'plan_sha256':sha(A/'plan.json'),'DE_probability_evaluation':'Average precision on FDR-test-eligible genes only, own coordinate excluded; untested genes are not negative labels','LFC_evaluation':'Conditional error/direction on actual DE genes only. Model probability is not a p-value.'})
 print(frame.groupby(['source','seen'])[['raw_pDE_AP_eligible_only','DE_prevalence','raw_PIE_LFC_DE_direction','V24_intended_DE_direction','calibrated_intended_DE_direction']].mean().to_string(),flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
