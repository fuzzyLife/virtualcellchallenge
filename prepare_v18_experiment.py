"""Create controls-only PIE inputs; freeze all comparisons before any outcome is scored."""
import json
import numpy as np
import build_submission as b
from v18_methods import ROOT,A,FAMILIES,NAMES,shuffled_indices
from fit_v7_response import raw_inputs,read_cache

def main(inputs_only=False):
 A.mkdir(parents=True,exist_ok=True)
 if (A/'plan.json').exists():raise FileExistsError('Preserve V18 freeze')
 data=read_cache();old=json.loads((ROOT/'experiments/v17_learned_heterogeneity/analysis/plan.json').read_text());confirm=json.loads((ROOT/'experiments/v17_learned_heterogeneity/analysis/confirmation_plan.json').read_text())
 train=json.loads((ROOT/'analysis/post_v17_research/PIE_training_meta.json').read_text());gt=json.loads((ROOT/'reference/v18_pie/sources/gene_text/meta.json').read_text());axis=[g for g in data['genes'] if g in set(train['genes'])];assert len(axis)==6290 and set(axis)<=set(gt['keys']);assert not set(data['targets'])&set(train['pert_to_id'])
 inputs=A/'inputs';inputs.mkdir(exist_ok=True);positions=np.array([list(data['genes']).index(g) for g in axis]);contexts=[]
 for ctx,ver in zip(old['contexts'],confirm['contexts']):
  assert ctx['family']==ver['family'];c=dict(ctx);c['context']=c['family']+'_v18';c['cells']=[r for r in ctx['cells'] if r['split']=='development']+[r for r in ver['cells'] if r['split']=='confirmation'];assert not {r['target'] for r in c['cells'] if r['split']=='development'}&{r['target'] for r in c['cells'] if r['split']=='confirmation'}
  x,meta,measured,source=raw_inputs(c['family'],data['genes']);ids=np.array(c['left']);assert np.all(meta.gene_target.to_numpy(str)[ids]=='Non-Targeting')
  depth=meta.total_counts.to_numpy(float)[ids];assert np.all(depth>=np.asarray(x[ids].sum(1)).ravel()*.999)
  y=x[ids].astype(np.float32).multiply((1e4/depth).astype(np.float32)[:,None]).tocsr();y.data=np.log1p(y.data);mean=np.asarray(y.mean(0)).ravel()[positions];assert np.isfinite(mean).all()
  file=inputs/(c['family']+'.npz');np.savez_compressed(file,genes=np.asarray(axis,str),competition_genes=data['genes'],targets=data['targets'],ctrl_mean=mean,measured=measured[positions],positions=positions,control_rows=ids)
  c['pie_controls_sha256']=b.sha256(file);contexts.append(c)
 if inputs_only:
  print('Prepared controls-only V18 inputs',flush=True);return
 deps=['collect_v18_pie_assets.py','v18_methods.py','prepare_v18_experiment.py','infer_v18_pie.py','run_v18_experiment.py','assess_v18_experiment.py','test_v18_pie.py','audit_v18_experiment.py','v16_methods.py','fit_v7_response.py','v12_cpu_metrics.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py']
 record={'version':18,'contexts':contexts,'seeds':old['seeds'],'candidates':NAMES,'eligible':['pie25','confidence'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'assets_sha256':b.sha256(A/'assets.json'),'parent_plan_sha256':b.sha256(ROOT/'experiments/v17_learned_heterogeneity/analysis/plan.json'),'confirmation_parent_sha256':b.sha256(ROOT/'experiments/v17_learned_heterogeneity/analysis/confirmation_plan.json'),'source_cache_sha256':old['source_cache_sha256'],'training_metadata_sha256':b.sha256(ROOT/'analysis/post_v17_research/PIE_training_meta.json'),'shuffle':shuffled_indices(300).tolist(),
 'expert':'PIE K562-held-out checkpoint. Strict full state_dict load, FP32, batch1; 6290 measured training-vocabulary outputs only. All 300 targets absent from training. Unknown context text omitted. Original per-source features; every original NCBI token for targets retained.',
 'evidence':'Novel query dataset names not present in replogle training and targets absent from entire training pert_to_id: all seven evidence blocks exactly zero by author serve_evidence semantics. Preserve original split hash/evidence key; unit-test missing-record equivalence rather than download irrelevant response cache.',
 'controls':'Mean of natural log1p per-cell CP10k; original full source library sizes before gene-axis restriction. Only left NTC used. No treatment values in inference.',
 'blend':'PIE log2 FC clipped -3..3 then same V16 expression gate. pie25 weights .25 on covered measured coordinates, confidence weights .25*clip((p_de-.5)/.5,0,1). p_de is not FDR. Own coordinate unchanged; V16 residual10 emitter and all seeds fixed. amplitude matches confidence full-axis non-own L1; shuffled uses fixed derangement of PIE rows before same gating/blending.',
 'release':'Existing raw/worst-family guards; >=3 available metrics improve including DE; confidence must beat shuffled and amplitude on >=2 of fidelity/reach/Jaccard. pie25 likewise must beat both controls. Missing NMAE remains missing. Freeze once-only verification on same rule, never retune after failure. No ranking guarantee.',
 'verification':'Reuse 8 corrected, still-unopened V17 target/family pairs (2 seeds). Algorithmically reserved, not blind biology; older public outcomes inspected in earlier work. Do not score unless development passes.'}
 b.write_json(A/'plan.json',record);print('V18 frozen',b.sha256(A/'plan.json'),flush=True)
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--inputs-only',action='store_true');main(parser.parse_args().inputs_only)
