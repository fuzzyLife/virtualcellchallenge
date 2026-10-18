"""Freeze V57 before evaluation; build control-only queries and disjoint-view labels."""
import argparse,gc
import numpy as np
from v57_pie import *

def freeze():
 import pandas as pd
 assert not (A/'plan.json').exists(),'Preserve the V57 plan'
 targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.astype(str).tolist()
 seen=set(read(ROOT/'analysis/v56_pie_release/audit.json')['training_split_coverage']['xd']['competition_targets'])
 perm=np.arange(300)
 for value in [True,False]:
  ids=np.array([i for i,t in enumerate(targets) if (t in seen)==value]);order=np.random.default_rng(SEED+int(value)).permutation(ids);perm[order]=np.roll(order,1)
 assert np.all(perm!=np.arange(300))
 cohorts={};inputs={}
 for source in SOURCES:
  f=ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(source+'.json');d=read(f);inputs[str(f.relative_to(ROOT))]=sha(f)
  groups=d['groups'];selected=[]
  if True:
   for value in [True,False]:
    group=[g for g in groups if (g['target'] in seen)==value]
    group.sort(key=lambda g:hashlib.sha256(('v57-pilot:'+source+':'+g['target']).encode()).hexdigest())
    selected.extend(g['target'] for g in group[:(4 if source.startswith('CD4') else 2)])
  cohorts[source]={'targets':[g['target'] for g in groups],'seen_count':sum(g['target'] in seen for g in groups),'count_pilot':selected,'calibration_source':not source.startswith('CD4')}
 deps=['v57_pie.py','collect_v57_pie_assets.py','build_v57_pie_evidence.py','prepare_v57_pie.py','infer_v57_pie.py','fit_v57_pie_decoder.py','evaluate_v57_pie_counts.py','audit_v57_pie.py']
 cfg={'seed':SEED,'ridge_alphas':[10.,100.,1000.],'gene_samples_per_training_target':768,'mean_residual_bound':1.,'rank_bound':.25,'decoder_rank_floor':.005,'count_seeds':[57057,157057], 'arms':['v24','pie_mean','pie_joint','shuffled_joint'],'calibration_sources':['H1','HCT116','HEK293T'],'transfer_sources':['CD4_Rest_D2','CD4_Stim48_D2'],'maximum_real_treated_cells':400,'prediction_cells':400,'maximum_reference_controls':2000,'donor_split_sha256':'8d050f28aebf482135d42f758f6b440e0333b3dec8efc0dddab9dcc3099d4be4','GPU_precision':'FP32','inference_chunk':512}
 record={'version':57,'cohorts':cohorts,'seen_targets':sorted(seen),'unseen_targets':sorted(set(targets)-seen),'shuffle':perm.tolist(),'config':cfg,'inputs':inputs,'dependencies':{n:sha(ROOT/n) for n in deps},'leakage_contract':'No CD4 treatment labels fit or select decoder parameters. H1/Orion fitting is source-overlap calibration, not unseen-study evidence. Matching source target responses are subtracted from donor evidence. D3 remains unopened.','calibration_selection':'Whole-target split hash mod5==0 selects ridge alpha on H1/Orion only. Refit chosen alpha on all H1/Orion targets. Same fit strategy for trained shuffled-head control and evidence-only ablation. No selection from CD4 outcomes.','pilot_selection':'Four metadata-eligible targets per seen/unseen class per CD4 context, chosen by frozen target hashes, without scoring results. All available distinct actual treatment cells up to400; do not duplicate to inflate sample size.','decoder':'Positive count mixture of compact Poisson and empirical control residuals, with gene-wise residual weight from predicted rank effect and predicted mean change under a Gaussian approximation. Global dual moments and row library constraints refit before actual DE evaluation. Approximation is tested, not asserted exact.','release':'Experimental comparison only. No new submission unless held-CD4 mean/profile/DE evidence supports it. Report failures and missing metrics; do not replace them with zeros.'}
 write(A/'plan.json',record);print('V57 plan frozen',sha(A/'plan.json'),flush=True)

def prepare(source):
 import v52_joint_response as old
 import v38_response_hurdle as hist
 from v25_methods import raw,baseline
 from scipy import sparse
 p=plan();out=CACHE/'inputs';out.mkdir(exist_ok=True);path=out/(source+'.npz')
 if path.exists():assert sha(path)==read(A/'inputs'/(source+'.json'))['sha256'];return
 genes=old.data()['genes'];targets=old.data()['targets'];z=old.labels(source);x,meta,measured,receipt=raw(source,genes)
 labels=meta.gene_target.to_numpy(str);ntc=np.flatnonzero(labels=='Non-Targeting');assert len(ntc)>400
 summed=x.sum(1).A1.astype(float)
 if 'total_counts' in meta and np.isfinite(meta.total_counts.to_numpy(float)).all() and np.all(meta.total_counts.to_numpy(float)>=summed*.999):depth=meta.total_counts.to_numpy(float);depth_kind='Original metadata full-source library'
 else:depth=summed;depth_kind='All available aligned measured genes; original full-source library unavailable'
 def logcp(rows):
  y=x[np.asarray(rows)].astype(np.float64).multiply((1e4/depth[np.asarray(rows)])[:,None]).tocsr();y.data=np.log1p(y.data);return y
 ctrl=logcp(ntc).mean(0).A1
 group='CD4' if source.startswith('CD4') else ('Orion' if source in ['HCT116','HEK293T'] else 'H1')
 base=baseline(old.data(),z['query'],hist.legacy((group,)),production=True)
 truth=np.full((300,2,len(genes)),np.nan,np.float32)
 groups=read(ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(source+'.json'))['groups']
 for g in groups:
  ti=g['target_index'];n=np.array(g['cells_per_view'])
  for lane in g['lanes']:
   for view in range(2):
    tr=lane['treated_views'][view];cr=lane['control_views'][view]
    assert tr and cr and not set(tr)&set(cr) and np.all(labels[tr]==g['target']) and np.all(labels[cr]=='Non-Targeting')
    delta=(logcp(tr).mean(0).A1-logcp(cr).mean(0).A1)*len(tr)/n[view]
    if np.isnan(truth[ti,view]).all():truth[ti,view]=0
    truth[ti,view]+=delta.astype(np.float32)
 gt=set(read(R/'sources/gene_text/meta.json')['keys']);union=set()
 for q in (CACHE/'data').glob('*/meta.json'):union.update(read(q)['genes'])
 support=measured&np.isin(genes,list(gt&union));assert support.sum()>6000
 np.savez_compressed(path,genes=genes,targets=targets,ctrl_mean=ctrl.astype(np.float32),support=support,measured=measured,baseline=base,query_cpm=z['query'],mean=z['mean'],rank=z['rank'],valid=z['valid'],delta=truth)
 write(A/'inputs'/(source+'.json'),{'sha256':sha(path),'plan_sha256':sha(A/'plan.json'),'source':source,'source_record':receipt,'depth_definition':depth_kind,'NTC_rows':ntc.tolist(),'supported_genes':int(support.sum()),'treatment_views':'Exact V52 disjoint raw-row views, retained unchanged; input control expression uses NTC only','baseline':'V24 public hold-study implementation; recipient whole study excluded from its source fit; same baseline as V52 held-source comparisons','D3_opened':False})
 print('V57 query prepared',source,'genes',support.sum(),'targets',len(groups),flush=True);del x;gc.collect()

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['freeze','prepare']);q.add_argument('--source',choices=SOURCES);z=q.parse_args()
 if z.stage=='freeze':freeze()
 else:
  for s in ([z.source] if z.source else SOURCES):prepare(s)
