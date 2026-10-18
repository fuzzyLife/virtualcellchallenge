"""V60: V59 decomposition with NTC-only nuisance-adjusted state coordinates.

Depth/detection may contain biology; this is an explicit alternative tested
against the original decomposition, not an identified technical correction.
"""
from functools import lru_cache
import argparse
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v59_state_composition as m
from run_v59_transfer_stage import activate as activate_parent
PARENT_A=m.A;PARENT_CACHE=m.CACHE
A=m.ROOT/'experiments/v60_depth_adjusted_states/analysis'
CACHE=m.Path(m.tempfile.gettempdir())/'vcc_v60_depth_adjusted_states'


def residualize(values,ntc,cov):
 center=cov[ntc].mean(0);scale=np.maximum(cov[ntc].std(0),.1);z=np.clip((cov-center)/scale,-5,5)
 design=np.column_stack([np.ones(len(z)),z,z*z,z[:,0]*z[:,1]])
 c=design[ntc];penalty=np.eye(c.shape[1])*.1;penalty[0,0]=0
 coef=np.linalg.solve(c.T@c+penalty,c.T@values[ntc])
 residual=values-design@coef;sd=np.maximum(residual[ntc].std(0),.5)
 result=np.clip(residual/sd,-5,5).astype(np.float32)
 return result,dict(covariate_center=center.tolist(),covariate_scale=scale.tolist(),coef=coef.tolist(),residual_scale=sd.tolist())


def standardize(x,meta,positions):
 # Neither nuisance covariate contains any of the 300 direct-target coordinates.
 nonpanel=~np.isin(m.axes()[0],m.axes()[1]);np_counts=x[:,nonpanel]
 library=np.maximum(np_counts.sum(1).A1/1e6*meta['depth'],1.)
 detected=np_counts.getnnz(1);cov=np.column_stack([np.log(library),np.log1p(detected)])
 result=np.zeros((x.shape[0],len(positions)),np.float32);traces=[]
 for lane in sorted(set(meta['lane'])):
  rows=np.flatnonzero(meta['lane']==lane);ntc=np.flatnonzero(meta['label'][rows]=='Non-Targeting');assert len(ntc)>=2
  values=np.log1p(x[rows][:,positions].toarray());v,trace=residualize(values,ntc,cov[rows]);result[rows]=v
  trace.update(lane=str(lane),control_rows=meta['rows'][rows[ntc]].tolist());traces.append(trace)
 return result,traces


def freeze():
 assert not (A/'plan.json').exists()
 p=m.read(PARENT_A/'plan.json');p['version']=60;p['parent_plan_sha256']=m.sha(PARENT_A/'plan.json')
 p['parent_CPM_correction_sha256']=m.sha(PARENT_A/'pre_score_correction.json');p['parent_mask_correction_sha256']=m.sha(PARENT_A/'pre_transfer_mask_correction.json')
 p['dependencies']['v60_depth_adjusted_states.py']=m.sha(__file__);p['dependencies']['run_v60_experiment.py']=m.sha(m.ROOT/'run_v60_experiment.py');p['dependencies']['test_v60_states.py']=m.sha(m.ROOT/'test_v60_states.py')
 p['state_adjustment']='NTC-only lane-specific ridge of log1p(CPM) on centered/scaled log nonpanel raw depth, log nonpanel detected-gene count, their squares and interaction; intercept unpenalized, other ridge .1. Divide residual by NTC residual SD floor .5 and clip +/-5. All 300 panel genes excluded from nuisance covariates. Source fold gene selection, PCA, KMeans, response/shuffle/decoder methods and release gates unchanged.'
 p['reason']='V59 NTC-only quality audit finds CD4 PC1/PC2 correlated .87-.92 with library and .89-.94 with detection. This independent version tests whether nuisance-driven alignment obscures transferable state mechanisms. It may also remove biological cell-size/complexity effects; positive results must survive end-to-end tests.'
 p['version_selection']='V60 primary joint50 and gates fixed before V60 state/treatment outcomes. Compare with V59 descriptively; do not select diagnostic arms post hoc. No claim that all original abundance structure was technical.'
 p['raw_cache_reuse']={s:m.sha(PARENT_A/'raw'/(s+'.json')) for s in m.SOURCES};m.write(A/'plan.json',p);print('V60 frozen',m.sha(A/'plan.json'),flush=True)


def activate():
 activate_parent();m.A=A;m.CACHE=CACHE
 @lru_cache(1)
 def plan():
  p=m.read(A/'plan.json')
  for n,h in p['dependencies'].items():assert m.sha(m.ROOT/n)==h,'V60 dependency changed: '+n
  return p
 @lru_cache(5)
 def data(source):
  p=plan();rpath=PARENT_A/'raw'/(source+'.json');assert m.sha(rpath)==p['raw_cache_reuse'][source];r=m.read(rpath);path=PARENT_CACHE/'raw'/(source+'.npz');side=path.with_name(source+'_meta.npz');assert m.sha(path)==r['sha256'] and m.sha(side)==r['meta_sha256'];return sparse.load_npz(path),m.load(side)
 m.plan=plan;m.source_data=data;m.standardize=standardize
 # Parent closure wrapper is already installed. Its model reader now resolves V60 paths;
 # measured-mask routing stays installed as well. No V59 file is rewritten.
 return m

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','states','mechanisms','transfer','evaluate','count','count_summary','plot']);p.add_argument('--source',choices=m.SOURCES);p.add_argument('--outer',choices=(*m.OUTERS,'all'));a=p.parse_args()
 if a.stage=='freeze':freeze()
 else:
  activate()
  with threadpool_limits(limits=2):
   if a.stage=='states':m.fit_states(a.outer)
   elif a.stage=='mechanisms':m.mechanisms(a.outer,a.source)
   elif a.stage=='transfer':m.transfer(a.source)
   elif a.stage=='evaluate':m.evaluate()
   elif a.stage in ['count','count_summary']:
    import evaluate_v59_counts as counts
    if a.stage=='count':counts.run(a.source)
    else:counts.summarize()
   else:
    import plot_v59_states as plot
    # The original plot module imports its parent's activate symbol. Bind the
    # V60 activation instead so its source methods and output directory agree.
    plot.activate=lambda:None;plot.main()
