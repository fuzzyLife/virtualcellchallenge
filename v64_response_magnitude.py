"""V64: calibrate downstream response magnitude, retaining V24 gene directions."""
from pathlib import Path
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v62_conserved_response as parent
from v57_pie import read,write,sha
import build_submission as b
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v64_response_magnitude';A=OUT/'analysis';ARMS=('magnitude',)
data=parent.data

def calibrate(base,teacher,genes,targets,shuffle=False):
 axis=~np.isin(genes,targets);n0=np.linalg.norm(base[:,axis],axis=1);n1=np.linalg.norm(teacher[:,axis],axis=1)
 scale=np.clip(np.divide(n1,n0,out=np.ones(len(n0)),where=n0>0),0,2)
 if shuffle:scale=np.roll(scale,1)
 result=np.clip(base*scale[:,None],-3,3).astype(np.float32);own=np.array([list(genes).index(t) for t in targets]);result[np.arange(len(targets)),own]=0
 assert np.all(result[base==0]==0) and np.all(result*base>=0)
 return result,scale

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v64_response_magnitude.py','run_v64_counts.py','test_v64_magnitude.py','run_v62_counts.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 write(A/'plan.json',dict(candidate='magnitude',dependencies={n:sha(ROOT/n) for n in deps},parent_plan_sha256=sha(parent.A/'plan.json'),teacher_plan_sha256=sha(ROOT/'experiments/v63_identity_guard/analysis/plan.json'),known_outcomes='V62/V63 means improve but count discrimination fails. A post-development 44-case desired-mean diagnostic finds V24 LFC norm calibration improves both mean error and PDS. V64 separately tests final integer cells; this is outcome-informed, repeatedly reused public development. No previous plan is amended.',estimator='Scale the V24 downstream LFC vector by guarded-teacher/V24 nonpanel LFC norm ratio per target/context, clipped0..2. Preserve all signs/zeros, clip LFC+/-3, keep own LFC zero and existing emitter own factor0.1. No transfer of teacher gene locations/signs. Norm is not observed knockdown efficacy and is not exactly the nonlinear assay-space norm.',negative_control='Cyclically derange the 300 inferred scalar magnitudes within the same context; retain V24 gene-response identity in both arms. Tests usefulness of target-specific magnitude allocation, not the entire V24 predictor.',selection='Only fixed magnitude primary. Same20 development count cases; retain prior six-metric/PDS/FID/REACH/JAC/NMAE/expression guards and beats-shuffle expression/PDS requirements. Remaining24 predefined cases at two seeds only after development support; otherwise label any extra reused-public evaluation diagnostic. No candidate/parameter reselection. No private anchors or new-study claim.',D3_opened=False,automatic_upload=False))
 print('V64 frozen',sha(A/'plan.json'),flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert sha(parent.A/'plan.json')==p['parent_plan_sha256'];return p

def prepare():
 check();from run_v35_v36 import contexts
 from run_v37_ensemble import checked_folder
 d=data();dest=A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family'];path=dest/(name+'.npz');source=ROOT/'experiments/v63_identity_guard/analysis/predictions'/(name+'.npz');r=read(source.with_suffix('.json'));assert sha(source)==r['sha256']
  with np.load(source,allow_pickle=False) as z:base=z['all_baseline'].copy();teacher=z['all_effects'].copy();query=z['query'].copy();measured=z['measured'].copy()
  actual,scale=calibrate(base,teacher,d['genes'],d['targets']);shuffled,_=calibrate(base,teacher,d['genes'],d['targets'],True);ix=np.array([list(d['targets']).index(c['target']) for c in ctx['cells']])
  np.savez_compressed(path,query=query,measured=measured,v24=base[ix],magnitude=actual[ix],magnitude_shuffled=shuffled[ix],all_effects=actual,all_baseline=base,scales=scale)
  write(path.with_suffix('.json'),dict(sha256=sha(path),teacher_prediction_sha256=r['sha256'],minimum_scale=float(scale.min()),maximum_scale=float(scale.max()),median_scale=float(np.median(scale)),all_signs_and_zeros_preserved=True));print('V64 prepared',name,'median factor',float(np.median(scale)),flush=True)
 write(A/'mean_results.json',dict(selected='magnitude',selection='Only fixed primary from separately frozen V64 plan; no count-outcome selection',teacher_mean_control_sha256=sha(ROOT/'experiments/v63_identity_guard/analysis/amplitude_control_results.json')))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);args=p.parse_args()
 with threadpool_limits(limits=2):globals()[args.stage]()
