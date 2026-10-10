"""Independently audit V12 cell identities, split boundaries, source hashes, and repeated seeds."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def main(folder):
 path=A/folder/'validation_plan.json';plan=json.loads(path.read_text());assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256']
 assignments={'development':set(),'confirmation':set()};rows=[];sources={}
 for family in ['H1','HCT116','HEK293T']:
  directory=ROOT/('reference/v5_h1' if family=='H1' else f'reference/v7_sources/{family}');file=directory/'selected_cells.csv';rec=json.loads((directory/'extraction.json').read_text());expected=rec['selected_cells_sha256' if family=='H1' else 'selected_sha256'];assert build.sha256(file)==expected
  meta=pd.read_csv(file)
  if family=='H1':meta['sample']=meta.batch.astype(str);meta['gene_target']=meta.target_gene.replace({'non-targeting':'Non-Targeting'});meta['cell_id']=meta.source_cell_id
  meta['sample']=meta['sample'].astype(str);assert meta.cell_id.is_unique;sources[family]={'metadata_sha256':expected,'cells':len(meta)}
  for ctx in [c for c in plan['contexts'] if c['family']==family]:
   left,right=set(ctx['left']),set(ctx['right']);assert not left&right
   assert meta.iloc[list(left|right)].gene_target.eq('Non-Targeting').all()
   assert meta.iloc[list(left|right)]['sample'].isin(ctx['batches']).all()
   for seed in plan['seeds']:
    cells=[c for c in ctx['cells'] if c['seed']==seed]
    assert all(c['template']==cells[0]['template'] and c['reference']==cells[0]['reference'] for c in cells)
    for c in cells:
     template,reference,actual=map(set,[c['template'],c['reference'],c['actual']]);assert len(template)==len(c['template'])==400 and len(reference)==len(c['reference'])==400
     assert template<=left and reference<=right and not actual&(left|right)
     assert 2<=len(actual)==len(c['actual'])<=400 and meta.iloc[c['actual']].gene_target.eq(c['target']).all() and meta.iloc[c['actual']]['sample'].isin(ctx['batches']).all()
     assignments[c['split']].add(c['target']);rows.append({'family':family,'context':ctx['context'],'target':c['target'],'split':c['split'],'seed':seed,'actual_cells':len(actual),'template_cells':len(template),'reference_cells':len(reference),'full_reference_cells':len(right)})
 assert not assignments['development']&assignments['confirmation']
 pd.DataFrame(rows).to_csv(A/folder/'split_identity_cases.csv',index=False)
 build.write_json(A/folder/'split_identity_audit.json',{'passed':True,'case_count':len(rows),'sources':sources,'plan_sha256':build.sha256(path),'script_sha256':build.sha256(Path(__file__)),'source_cache_hash_verified':True,'globally_disjoint_target_sets':True,'target_independent_shared_templates':True,'actual_cells_unique_with_correct_target_batch':True,'template_reference_actual_identity_disjoint':True,'limitation':'Seeds reuse some or all actual treatment cells, and technical strata are not independent biological families. Do not treat seeds or individual cells as independent experiments.'})
 print('Split identity and input audit passed:',folder,len(rows),'cases',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--folder',required=True);main(p.parse_args().folder)
