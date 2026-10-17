"""Bind the trained V53 checkpoint, model API and completed scientific audits."""
from pathlib import Path
import numpy as np,torch
from threadpoolctl import threadpool_limits
import run_v53_joint_guard as guard
import v52_joint_response as m

def main():
 guard.install();m.check()
 names=['label_audit.json','boundary_audit.json','feature_contract_audit.json','final_audit.json']
 for name in names:assert m.read(m.A/name)['passed'],name
 from predict_v52_joint_response import predict
 z=m.labels('H1');ti=int(np.flatnonzero(z['valid'].any(1))[0]);t=z['targets'][ti]
 out=predict(z['query'],outer='all',targets=[t],device='cpu')
 assert out['mean'].shape==out['rank'].shape==out['valid'].shape==(1,18533)
 assert np.isfinite(out['mean']).all() and np.isfinite(out['rank']).all()
 assert abs(out['mean']).max()<=3 and abs(out['rank']).max()<=.5
 own=z['genes']==t;assert not out['valid'][0,own].any() and (out['mean'][0,own]==0).all() and (out['rank'][0,own]==0).all()
 model=m.A/'models/all/joint.pt';rec=m.read(model.with_suffix('.json'));assert m.b.sha256(model)==rec['sha256']
 evaluation=m.read(m.A/'evaluation.json')
 proof=dict(passed=True,model_kind='Experimental all-source joint mean/rank response predictor',checkpoint=str(model.relative_to(m.ROOT)),checkpoint_sha256=rec['sha256'],best_source_validation_epoch=rec['best_epoch'],
  api_file='predict_v52_joint_response.py',api_sha256=m.b.sha256(m.ROOT/'predict_v52_joint_response.py'),
  install_file='run_v53_joint_guard.py',install_sha256=m.b.sha256(m.ROOT/'run_v53_joint_guard.py'),
  guarded_plan_sha256=m.b.sha256(m.A/'plan.json'),parent_label_plan_sha256=m.b.sha256(guard.BASE/'plan.json'),
  audit_receipts_sha256={name:m.b.sha256(m.A/name) for name in names},evaluation_sha256=m.b.sha256(m.A/'evaluation.json'),
  pilot_checks=evaluation['checks'],response_pilot_passed=evaluation['pilot_success'],upload_ready=False,D3_opened=False,
  smoke_test=dict(all_source_checkpoint_reloaded=True,device='cpu',query_source='H1 source controls (numerical smoke only, not a held-study test)',targets=1,genes=18533,finite_outputs=True,bounds_passed=True,own_gene_excluded=True),
  finalizer_sha256=m.b.sha256(Path(__file__)),scope='Audited trained response model; no count decoder, upload package, new blind test or private-score gain is asserted.')
 m.write(m.A/'delivery_manifest.json',proof);print('V53 all-source model/API smoke and delivery binding passed',flush=True)
if __name__=='__main__':
 torch.set_num_threads(2)
 with threadpool_limits(limits=2):main()
