"""Freeze the vocabulary-compatibility follow-up before its active joint tests."""
import json
import build_submission as b
from v21_methods import ROOT,A,NAMES,FAMILIES
import v20_methods as parent

def main():
 A.mkdir(parents=True,exist_ok=True)
 if (A/'plan.json').exists():raise FileExistsError('Preserve V21 freeze')
 old=parent.load_plan();deps=['v21_methods.py','fit_v21_programs.py','prepare_v21_experiment.py','infer_v21_stack.py','run_v21_experiment.py','assess_v21_experiment.py','audit_v21_experiment.py','test_v21_counts.py','v20_methods.py','stack_v20_mean_adapter.py','v19_methods.py','v16_methods.py','fit_v7_response.py','v12_cpu_metrics.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py'];coverage=[]
 for family in FAMILIES:
  rec=json.loads((parent.A/'experts'/f'{family}_development.json').read_text())
  for row in rec['records']:
   for pair in row['pairs']:coverage.append({'family':family,'target':row['target'],'kind':pair['kind'],'coverage':pair['program_coverage']})
 assert all(r['coverage']<.8 for r in coverage)
 record={'version':21,'contexts':old['contexts'],'seeds':old['seeds'],'family_order':FAMILIES,'candidates':NAMES,'eligible':['stack_joint'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(parent.A/'plan.json'),'parent_assets_manifest_sha256':b.sha256(parent.A/'assets.json'),'old_program_coverage':coverage,
 'reason':'All V20 joint residuals abstained because inherited program vocabulary coverage was below its frozen .8 rule. V21 fixes input compatibility by fitting the same16-program/top2048 control-only SVD on Stack-covered measured genes, excluding all300 panel genes. No threshold is relaxed; no mean-response model or emitter coefficient is retuned. Some V20 development mean outcomes were already visible; new joint outcomes are not. This is not a fresh blind development sample.',
 'mean':'Preserve all V20 lfc,shrink,raw,supported,shuffled_lfc,shuffled_shrink response arrays bitwise. Recompute the exact same conditional mean outputs and check original generated-mean array hashes before projecting into compatible programs.',
 'joint':'Same source inputs, source size weights, conditional residual bounds, .10 empirical component, nearest query-state mapping, own-target handling and full integer moment constraints as V20. Add an architecture-matched shuffled-joint control with shuffled mean and shuffled conditional programs. Reweighting uses the same new programs.',
 'release':'Same available raw/worst-family guards and >=3 improved metrics including DE. Joint must beat amplitude, matched shuffled-joint, unchanged Stack mean and empirical reweighting on >=2 fidelity/reach/Jaccard; centered shape Wasserstein must also beat mean and reweighting. Only after passing may unspent8 reserved pairs be tested once. If V21 fails, retain a V20 mean finalist only if that separately frozen experiment passed. No private score promise.',
 'pretraining_overlap':old['pretraining_overlap'],'verification':old['verification']}
 b.write_json(A/'plan.json',record);print('V21 compatibility experiment frozen',b.sha256(A/'plan.json'),flush=True)
if __name__=='__main__':main()
