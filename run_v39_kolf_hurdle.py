"""Run V39 sequentially with the frozen V38 recipe and an added data source."""
import argparse,importlib.metadata,platform
from threadpoolctl import threadpool_limits
import v39_kolf_hurdle as extra
import v38_response_hurdle as m
import run_v38_response_hurdle as r
from v38_binning_adapter import install as install_binning
def freeze():
 extra.A.mkdir(parents=True,exist_ok=True);assert not (extra.A/'plan.json').exists();parent=m.read(m.ROOT/'experiments/v38_response_hurdle/analysis/plan.json')
 inputs=dict(parent['inputs']);inputs.update({str(p.relative_to(m.ROOT)):m.b.sha256(p) for p in [extra.A/'collection.json',extra.A/'response_features.npz',extra.A/'source_agreement.json',m.ROOT/'experiments/v38_response_hurdle/analysis/mean_results.json',m.ROOT/'experiments/v38_response_hurdle/analysis/specificity/results.json']})
 deps=dict(parent['dependencies']);deps.update({s:m.b.sha256(m.ROOT/s) for s in ['v39_kolf_hurdle.py','run_v39_kolf_hurdle.py','collect_v39_kolf_features.py','v38_binning_adapter.py']})
 m.write(extra.A/'plan.json',{**parent,'dependencies':deps,'inputs':inputs,'parent_V38_plan_sha256':m.b.sha256(m.ROOT/'experiments/v38_response_hurdle/analysis/plan.json'),'additional_hypothesis':'Seven auxiliary KOLF log-expression NTC-z features: signed/absolute/clipped value, product/agreement with pooled LFC, availability and quantization saturation. No label, dose, fold, target, gene sample, recipe, correction coefficient, generator or threshold change from V38.','data_scope':'Pinned all-perturbation browser export, 273 matching targets and18106 genes; not strong-only cells. No raw cells, variance or direct LFC conversion. KOLF is a separate public study in every source fold; H1 and KOLF share stem-cell biology, so whole-study exclusions are not independent-lineage guarantees.','known_outcomes':'V38 source/mean and zero-head ablation outcomes known, V38 raw treatment scoring in progress; original V38 immutable. Public diagnostics reused, no blind claim.','seed':m.SEED})
 print('V39 added-source plan frozen',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train','prepare','means','null','development','select','confirmation']);p.add_argument('--outer',choices=list(m.OUTERS)+['all']);a=p.parse_args();extra.install();install_binning()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  elif a.stage=='train':
   r.check()
   for o in [a.outer] if a.outer else m.OUTERS:m.fit(o)
  elif a.stage in ['null','development','confirmation']:r.counts(a.stage)
  else:getattr(r,a.stage)()

