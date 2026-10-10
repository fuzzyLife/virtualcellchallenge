"""Seal a completed V13 release without uploading or changing its predictor."""
import json,py_compile
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import build_submission as build
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v13_response_programs';A=OUT/'analysis'
V11_SHA='7841429ef987c43abcfb01a5be8f09864c5f8955e1bd7aae68fb4a03240be31d'

def read(name):return json.loads((A/name).read_text(encoding='utf-8'))

def main():
 ready=read('submission_ready.json');verified=read('independent_validation.json');model=read('reference_model.json')
 assert verified['passed'] and ready['official_prep_passed'] and not ready['uploaded']
 assert ready['sha256']==verified['package_sha256']==build.sha256(OUT/'prediction.vcc')
 assert model['confirmed_improvement']==read('confirmation_result.json')['confirmed_improvement']
 for name,digest in model['dependencies'].items():assert build.sha256(ROOT/name)==digest,name
 for name in ['response_checks.json','robust_source_checks.json','quality_source_checks.json','baseline_replay_audit.json']:
  record=read(name);assert record['passed'],name
 for folder in [A,A/'robust_stage',A/'quality_stage']:
  assert json.loads((folder/'split_identity_audit.json').read_text(encoding='utf-8'))['passed']
 old=ROOT/'experiments/v11_consensus_coverage/prediction.vcc'
 assert build.sha256(old)==V11_SHA,'Preserve the original published V11 package'
 scripts=sorted(ROOT.glob('*v13*.py'))
 for path in scripts:py_compile.compile(str(path),doraise=True)
 files=scripts+[ROOT/n for n in ['V13_STRATEGY.md','V13_SUBMISSION.md','v13_requirements.txt']]
 files += [p for p in A.rglob('*') if p.is_file() and p.suffix in {'.json','.csv','.py'} and 'count_cache' not in p.parts
           and p.name not in {'release_manifest.json','prediction_control_rows.csv'}]
 stages=[]
 for folder in [A,A/'robust_stage',A/'quality_stage']:
  frame=pd.read_csv(folder/'development_scores.csv')
  stages.append({'folder':str(folder.relative_to(A)),'case_candidate_populations':len(frame)//2,
    'candidates':sorted(frame.candidate.unique().tolist()),'reference_modes':sorted(frame.reference_mode.unique().tolist())})
 confirmation=pd.read_csv(A/model['validation_folder']/'confirmation_scores.csv')
 manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'package_sha256':ready['sha256'],
  'package_bytes':ready['bytes'],'validation_status':model['status'],'candidate':model['candidate'],
  'uploaded':False,'development':stages,'confirmation_case_candidate_populations':len(confirmation)//2,
  'count_population_definition':'One case/candidate population; its two reference modes reuse the same prediction. Baselines repeat across stages. Seeds, technical strata, and cells are not independent biological replication.',
  'published_v11_package_preserved_sha256':V11_SHA,'frozen_model_dependencies_verified':True,
  'independent_payload_audit_sha256':build.sha256(A/'independent_validation.json'),
  'full_sampling_provenance_sha256':build.sha256(A/'prediction_control_rows.csv'),
  'root_python_scripts_compiled':len(scripts),
  'tests':'Four meaningful tests across three response classes passed; exact held-family isolation and unknown-KD invariance. Their hashed outcome records are included.',
  'files':{str(p.relative_to(ROOT)).replace('\\','/'):build.sha256(p) for p in sorted(set(files))},
  'scope':'Integrity and reproducibility. Public confirmation and raw-guard status are separate from unknown leaderboard performance.'}
 build.write_json(A/'release_manifest.json',manifest)
 print('V13 release sealed:',model['status'],'package',ready['sha256'],flush=True)
if __name__=='__main__':main()
