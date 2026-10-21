"""Preserve existing scientific releases while testing V85/V86."""
from pathlib import Path
from v57_pie import read,write,sha
from seal_v69_release import paths
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v86_joint_posterior/analysis'
manifest=read(ROOT/'experiments/v69_count_bridge/analysis/release_manifest.json');expected={r['path'] for r in manifest['files']};actual={p.relative_to(ROOT).as_posix() for p in paths()};assert actual==expected
checked=0
for r in manifest['files']:
 p=ROOT/r['path'];assert p.stat().st_size==r['bytes']
 if p.suffix.lower() in ['.py','.json','.md','.csv','.txt']:assert sha(p)==r['sha256'];checked+=1
protocol=read(ROOT/'experiments/v84_intervention_operator/analysis/execution_protocol.json')
for name,h in {**protocol['code_sha256'],**protocol['inputs_sha256']}.items():assert sha(ROOT/name)==h,name
original=read(ROOT/'experiments/v85_guide_contrast/analysis/engineering_history/log_stage_plan.json')
for name,h in original['code_sha256'].items():
 p=ROOT/'experiments/v85_guide_contrast/analysis/engineering_history'/('log_stage_'+name+'.txt')
 if not p.exists():p=ROOT/name
 assert sha(p)==h,name
write(A/'preservation_audit.json',dict(passed=True,V69_inventory_files=len(expected),V69_metadata_code_hashes_checked=checked,V69_large_binary_hashes_rechecked=False,V84_frozen_code_and_recorded_inputs_unchanged=True,V85_original_code_and_plan_preserved=True,auditor_sha256=sha(__file__)))
print('Preservation audit passed',len(expected),checked,flush=True)
