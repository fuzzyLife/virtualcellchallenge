"""Pin the rejected-candidate evidence and verify the retained V11 package is unchanged."""
from pathlib import Path
import json,hashlib,py_compile
import numpy as np
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def main():
 choice=json.loads((A/'final_development_choice.json').read_text());assert choice['choice']['candidate']=='v11' and not choice['passing_choices'] and not choice['confirmation_read']
 assert not list(A.rglob('confirmation_scores.csv'))
 for folder,evidence in choice['evidence'].items():
  for name,h in evidence.items():assert build.sha256(A/folder/name)==h
 for folder in ['flex_stage','harmonized_stage']:
  assert not json.loads((A/folder/'development_assessment.json').read_text())['passing_candidates']
  assert json.loads((A/folder/'development_evidence_audit.json').read_text())['passed']
  assert json.loads((A/folder/'source_checks.json').read_text())['passed']
 for name in ['metric_tests.json','official_cpu_parity.json','cross_stage_baseline_replay.json']:assert json.loads((A/name).read_text())['passed']
 populations=0;casecounts={}
 for folder in ['', 'intervention_stage','flex_stage','harmonized_stage','poisson_stage']:
  cache=A/folder/'count_cache';cases=list(cache.glob('*_development.json'));casecounts[folder or 'initial_distribution']=len(cases)
  for side in cases:
   r=json.loads(side.read_text());assert build.sha256(side.with_suffix('.npz'))==r['sha256'];modes={row['reference_mode'] for row in r['records']};assert modes=={'sampled400','full_independent'}
   names={row['candidate'] for row in r['records']};assert len(r['records'])==2*len(names);populations+=len(names)
 assert populations==952
 extraction=json.loads((A/'flex_extraction.json').read_text());source=ROOT/'reference/v12_resources/flex_source'
 for name,key in [('summaries.npz','summary_sha256'),('selected_counts.npz','counts_sha256'),('selected_cells.csv','metadata_sha256')]:assert build.sha256(source/name)==extraction[key]
 for name,key in [('extract_v12_flex.py','script_sha256'),('extract_v12_flex.R','r_script_sha256')]:assert build.sha256(ROOT/name)==extraction[key]
 assert extraction['unique_cells']==22733 and extraction['unique_controls']==7949
 prior=json.loads((ROOT/'analysis/post_v11_research/research_manifest.json').read_text());v11=ROOT/'experiments/v11_consensus_coverage/analysis'
 for name,h in prior['frozen_v11_dependencies_verified'].items():assert build.sha256(ROOT/name)==h
 for name,h in prior['frozen_v11_snapshots_verified'].items():assert build.sha256(v11/name)==h
 model=json.loads((v11/'reference_model.json').read_text())
 for name,h in model['dependencies'].items():assert build.sha256(ROOT/name)==h
 ready=json.loads((v11/'submission_ready.json').read_text());package=ROOT/'experiments/v11_consensus_coverage/prediction.vcc';assert package.stat().st_size==ready['bytes'] and build.sha256(package)==ready['sha256']
 print('V11 package and all frozen snapshots unchanged',flush=True)
 code={p.name:build.sha256(p) for p in sorted(ROOT.glob('*v12*.py'))}
 for p in sorted(ROOT.glob('*v12*.py')):py_compile.compile(str(p),doraise=True)
 code.update({p.name:build.sha256(p) for p in sorted(ROOT.glob('*v12*.R'))});code['v12_official_requirements.txt']=build.sha256(ROOT/'v12_official_requirements.txt')
 records={}
 for p in sorted(A.rglob('*')):
  if not p.is_file() or 'count_cache' in p.parts or p.suffix in ['.npz','.h5ad'] or p.name.startswith('neutral_genes_') or p.name=='research_final_manifest.json':continue
  records[str(p.relative_to(A)).replace('\\','/')]=build.sha256(p)
 result={'passed':True,'decision':'Retain V11; all new count/source alternatives fail fixed guards. No V12 package, upload, or confirmation count evaluation.','generated_test_populations':populations,'completed_case_caches':casecounts,'retained_v11_package_sha256':ready['sha256'],'retained_v11_package_bytes':ready['bytes'],'prior_v11_dependencies_verified':len(model['dependencies']),'prior_v11_snapshots_verified':len(prior['frozen_v11_snapshots_verified']),'matched_Flex_source_hashes_verified':True,'root_code_sha256':code,'small_evidence_sha256':records,'report_sha256':build.sha256(ROOT/'V12_RESEARCH.md'),'limitations':['Reused public development biology and multiple bounded hypothesis stages; no fully blind biological claim.','Three actual-count source families; technical strata and seeds are not independent biological replicates.','One sampled-reference HEK panel has no NMAE-eligible target; complete-panel guard not relaxed.','Public source anchors can be noisy and are not private leaderboard anchors.','Poisson stage terminated after complete H1 worst-family guard failure; partial later cases retained.']}
 build.write_json(A/'research_final_manifest.json',result);print('Final V12 evidence pinned:',populations,'generated populations; no replacement promoted',flush=True)
if __name__=='__main__':main()
