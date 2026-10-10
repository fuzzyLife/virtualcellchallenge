"""Control-only tests of candidate count distributions; no submission model selection."""
import gc,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from v6_methods import moments
from v12_poisson_methods import emit
from v12_cpu_metrics import de
from v11_methods import SETTINGS
from evaluate_atlas_transfer import load_atlas_model
ROOT=build.ROOT;OUT=ROOT/'experiments/v12_distribution_calibrated/analysis/poisson_stage/neutral';SEEDS=[2026100911,2026100927]
METHODS=['real_resample','v11','noise85','noise100']

def simulate(template,depth,p,b,var,method,seed,atlas,fit_inverse):
 if method=='real_resample':return template.copy(),{'moment_matching':False}
 return emit(template,np.arange(len(depth)),np.zeros(len(p)),None,seed,method,prepared=(depth,p,b,var,fit_inverse),atlas=atlas)


def main():
 OUT.mkdir(parents=True,exist_ok=True);planpath=OUT/'neutral_test_plan.json'
 plan={'methods':METHODS,'seeds':SEEDS,'contexts':list('ABC'),'prediction_cells':400,'reference_cells':[400,2000],
  'split':'Per-guide disjoint 50% fit, 25% templates, 25% references; control-only; no own-target knockdown.',
  'script_sha256':build.sha256(Path(__file__)),'v11_model_sha256':build.sha256(ROOT/'experiments/v11_consensus_coverage/analysis/reference_model.json'),
  'purpose':'Distribution and null diagnostics only; no parameter sweep or winner selection; no hidden response labels.',
  'gamma_note':'Independent per-gene overdispersion is a diagnostic, not a scVI model; correlated biological structure is absent.',
  'empirical_note':'Full empirical composition template preserves covariance with both moments constrained; not a CellOT reproduction.',
  'dependencies':{n:build.sha256(ROOT/n) for n in ['v12_poisson_methods.py','v12_cpu_metrics.py','v11_methods.py','v6_methods.py','evaluate_atlas_transfer.py','build_submission.py']}}
 if planpath.exists():assert json.loads(planpath.read_text())==plan
 else:build.write_json(planpath,plan)
 manifest,genes,targets=build.read_spec(ROOT);atlas,atlas_sha=load_atlas_model();records=[];gene_records=[];inputs=[]
 for ci,context in enumerate('ABC'):
  x,guides,library,*_=build.load_control(ROOT,context,genes,manifest);inputs.append({'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
  for seed in SEEDS:
   side=OUT/f'neutral_{context}_{seed}.json';gene_path=OUT/f'neutral_genes_{context}_{seed}.csv'
   if side.exists():
    saved=json.loads(side.read_text());assert saved['plan_sha256']==build.sha256(planpath) and saved['gene_sha256']==build.sha256(gene_path)
    records.extend(saved['records']);print('Verified prior diagnostic',context,seed,flush=True);continue
   rng=np.random.default_rng(seed+ci);fit=[];templates=[];references=[]
   for guide in sorted(set(guides)):
    rows=rng.permutation(np.flatnonzero(guides==guide));a=len(rows)//2;b=3*len(rows)//4
    fit.extend(rows[:a]);templates.extend(rows[a:b]);references.extend(rows[b:])
   assert not set(fit)&set(templates) and not set(fit)&set(references) and not set(templates)&set(references)
   depth,p,b,var=moments(x[fit]);ids=build.guide_stratified_sample(guides[templates],library[templates],400,rng)
   template=x[np.asarray(templates)[ids]];libs=np.asarray(template.sum(axis=1)).ravel().astype(np.int64)
   refs=x[rng.choice(references,2000,replace=False)];refvar=moments(refs)[3];eligible=moments(refs)[1]>5e-6
   group=[];pergene=[]
   for mi,method in enumerate(METHODS):
    count,diag=simulate(template,libs,p,b,var,method,seed+ci,atlas,float(np.mean(1/depth)))
    _,cp,cb,cv=moments(count);detect=np.asarray((count>0).mean(axis=0)).ravel();rdetect=np.asarray((refs>0).mean(axis=0)).ravel()
    ratio=np.divide(cv,refvar,out=np.full(len(p),np.nan),where=refvar>0)
    base={'context':context,'seed':seed,'method':method,'CPM_variance_ratio':float(cv.sum()/refvar.sum()),
     'median_gene_variance_ratio':float(np.nanmedian(ratio[eligible])),'detection_MAE':float(abs(detect-rdetect)[eligible].mean()),
     'mean_probability_L1':float(abs(cp-p).sum()),'bulk_probability_L1':float(abs(cb-b).sum()),
     'mean_detected_genes':float(count.getnnz(axis=1).mean()),'exact_original_depths':bool(np.array_equal(np.asarray(count.sum(axis=1)).ravel(),libs)),
     'positive_integer_counts':bool(np.all(count.data>0)&np.all(count.data==np.floor(count.data)))}
    for size in [400,2000]:
     test=de(count,refs[:size],genes);n=int(test['eligible'].sum());calls=int((test['eligible']&(test['q']<.05)).sum())
     group.append({**base,'reference_cells':size,'null_DE_calls':calls,'eligible_genes':n,'null_call_fraction':calls/max(n,1)})
    for j in np.flatnonzero(eligible):pergene.append({'context':context,'seed':seed,'method':method,'gene_name':genes[j],
     'fitted_mean_CPM':float(p[j]*1e6),
     'variance_ratio':float(ratio[j]),'detection_error':float(detect[j]-rdetect[j])})
    print('Neutral diagnostic',context,seed,method,'null DE',group[-1]['null_DE_calls'],flush=True)
   pd.DataFrame(pergene).to_csv(gene_path,index=False)
   build.write_json(side,{'plan_sha256':build.sha256(planpath),'gene_sha256':build.sha256(gene_path),'records':group});records.extend(group)
  del x;gc.collect()
 frame=pd.DataFrame(records);frame.to_csv(OUT/'neutral_scores.csv',index=False)
 summary=frame.groupby(['method','reference_cells']).agg(max_null_DE_calls=('null_DE_calls','max'),mean_null_DE_calls=('null_DE_calls','mean'),
  max_null_call_fraction=('null_call_fraction','max'),mean_variance_ratio=('CPM_variance_ratio','mean'),median_gene_variance_ratio=('median_gene_variance_ratio','mean'),mean_detection_MAE=('detection_MAE','mean')).reset_index()
 summary.to_csv(OUT/'neutral_summary.csv',index=False)
 build.write_json(OUT/'neutral_test_manifest.json',{'completed':True,'plan_sha256':build.sha256(planpath),'inputs':inputs,'atlas_solver_sha256':atlas_sha,
  'records':len(frame),'scope':'Null-distribution pilot only, no perturbation performance claim. Two seeds and all three controls; references independent of fit and template populations.'})
 print(summary.to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
 frame=pd.read_csv(OUT/'neutral_scores.csv');assert frame.exact_original_depths.all() and frame.positive_integer_counts.all()
 build.write_json(OUT.parent/'source_checks.json',{'passed':True,'response':'Unchanged frozen V11; emission-only experiment','neutral_manifest_sha256':build.sha256(OUT/'neutral_test_manifest.json'),'exact_depths_and_integer_counts':True,'scope':'Neutral diagnostics report all contexts, including worst null-call fraction; do not claim a biological null-safety pass.'})
