"""Freeze V13 count comparison before reading generated-count outcomes."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v13_methods import ROOT
A=ROOT/'experiments/v13_response_programs/analysis/robust_stage'
from fit_v7_response import read_cache,raw_inputs
from v12_cpu_metrics import de
import build_submission as build

def main():
 path=A/'validation_plan.json';assert not path.exists()
 screen=json.loads((A/'response_screen_result.json').read_text(encoding='utf-8'))
 source=ROOT/'experiments/v12_distribution_calibrated/analysis/validation_plan.json'
 plan=json.loads(source.read_text(encoding='utf-8'))
 reserved=json.loads((A/'response_screen_plan.json').read_text(encoding='utf-8'))['reserved_targets']
 # Retain original disjoint target sets; new seeds do not make new biology.
 plan.update(version=13,emitters=['v11']+screen['advance'],reserved_targets=reserved,
 source_methods_sha256=build.sha256(ROOT/'v13_robust_methods.py'),cpu_metrics_sha256=build.sha256(ROOT/'v12_cpu_metrics.py'),
 response_screen_sha256=build.sha256(A/'response_screen_result.json'),
 parent_plan_sha256=build.sha256(source),script_sha256=build.sha256(Path(__file__)),
 source_response='Training-only robust source consensus; unchanged V11 noise65 integer emitter.',
 selection_policy='Six-metric whole-family-balanced primary full-independent references; sampled400 secondary with explicit missing NMAE. Original raw guards applied on available secondary metrics; primary must be complete. Choose greatest mean relative expression-error reduction among candidates passing guards, then assess once on reserved count confirmation. Failure means exploratory package only, never a claimed gain.',
 expression_improve_required_both_stages=True)
 # Preserve established seeds and controls so V11 replays can be checked against V12 caches.
 # Missing low-power secondary NMAE is reported, never imputed or silently called complete.
 build.write_json(path,plan)
 data=read_cache();genes=data['genes'];rows=[];assigned={'development':set(),'confirmation':set()}
 for family in ['H1','HCT116','HEK293T']:
  x,meta,_,_=raw_inputs(family,genes)
  for ctx in [c for c in plan['contexts'] if c['family']==family]:
   left,right=set(ctx['left']),set(ctx['right']);assert not left&right
   for cell in ctx['cells']:
    assert len(set(cell['actual']))==len(cell['actual']) and not set(cell['actual'])&(left|right)
    assert len(set(cell['template']))==len(set(cell['reference']))==400
    assert set(cell['template'])<=left and set(cell['reference'])<=right
    assert meta.iloc[cell['actual']].gene_target.eq(cell['target']).all()
    assigned[cell['split']].add(cell['target'])
  del x
 assert not assigned['development']&assigned['confirmation']
 build.write_json(A/'split_identity_audit.json',{'passed':True,'plan_sha256':build.sha256(path),'globally_disjoint_target_sets':True,'all_actual_rows_unique_and_target_correct':True,'template_and_reference_controls_disjoint':True,'scope':'Original V12 public development cases are reused; count-confirmation outcomes remain unopened. Seeds and technical strata are not independent biology.'})
 print('V13 count plan frozen:',plan['emitters'],flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
