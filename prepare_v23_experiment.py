"""Freeze the mean-preserving extension after V22 rejection, before V23 scores."""
import json,shutil
from pathlib import Path
from v23_methods import ROOT,A,NAMES,write
import v22_methods as parent
import build_submission as b

def main():
 A.mkdir(parents=True,exist_ok=True);dest=A/'plan.json'
 if dest.exists():raise FileExistsError('Preserve V23 frozen plan')
 old=parent.load_plan();selection=json.loads((parent.A/'selection.json').read_text());assert not selection['new_candidate_passes'] and not selection['checks']['nmae'];diagnostic=json.loads((A/'mean_drift_diagnostic.json').read_text());assert not diagnostic['treatment_DE_metrics_computed'] and diagnostic['mean_CPM_probability_L1_after']<diagnostic['mean_CPM_probability_L1_before']/100
 folder=A/'models';folder.mkdir(exist_ok=True);models={}
 for f in parent.FAMILIES+['CD4','all']:
  for suffix in ['.npz','.json']:
   name=f+suffix;shutil.copyfile(parent.A/'models'/name,folder/name);models[name]=b.sha256(parent.A/'models'/name);assert b.sha256(folder/name)==models[name]
 deps=['v23_methods.py','v23_mean_repair.py','prepare_v23_experiment.py','run_v23_experiment.py','assess_v23_experiment.py','test_v23_mean_preservation.py','v22_methods.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py']
 plan={'version':23,'contexts':old['contexts'],'names':NAMES,'eligible_finalists':['mean_preserved'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'parent_models':models,'parent_plan_sha256':b.sha256(parent.A/'plan.json'),'parent_failed_selection_sha256':b.sha256(parent.A/'selection.json'),'mean_drift_diagnostic_sha256':b.sha256(A/'mean_drift_diagnostic.json'),
 'mechanism':'Frozen V22 conditional populations plus reference-V16 integer moment repair. First match exact pooled gene totals using within-row moves; then bounded 2x2 count cycles preserve rows and columns while reducing mean-CPM differences. Max20k cycles, max2%library/200k transferred molecules; each improves weighted absolute mean error. No observed or hidden treatment truth used in repair.',
 'development':'Same14 public target/family pairs, two seeds, 400 real and predicted cells, six arms =168 full populations. Reused public development, not fresh biological confirmation. Only mean_preserved is a release finalist; conditional is the rejected V22 diagnostic.',
 'reserve':'Metadata-only corrected CD4 D2 Rest source (D1 only2 powered targets; D2 has15). Eight alphabetical targets>=400 with lane-matched NTC halves. All CD4 excluded from both mean bank and factor training. Model/guards fixed before any reserve scores. Same study as legacy CD4 summaries; new raw-cell check, not independent study. Reserve identity/extraction/mapping amendment saved under V22 analysis and copied before confirmation.',
 'release':'Identical V22 available-metric guard thresholds and >=3 improvements including DE, centered response-projection shape improves; beats scalar_preserved, independent_preserved and shuffled_preserved on >=2 DE metrics each. Null exactV16; new means match V16, gene totals/libraries exact. Same once-only reserve guards. Missing NMAE stays missing. No private anchor score estimated.',
 'interpretation':'A separately frozen correction of V22 emitted-mean drift, not retroactive relaxation of its failed NMAE guard. Underlying transferred V16 means can still be biologically wrong. If rejected, retain champion and label any upload-ready fallback honestly.'}
 write(dest,plan);print('Frozen V23 experiment',b.sha256(dest),flush=True)
if __name__=='__main__':main()
