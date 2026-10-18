"""Independent V61 population-label replay and empty lane-state accounting."""
import numpy as np
from scipy import sparse
import v61_pooled_state_response as v

def main():
 v.activate();m=v.m;genes,targets=m.axes();checks=[];empty=[]
 for source in m.SOURCES:
  x,meta=m.source_data(source);outer=m.GROUP[source];path=m.CACHE/'mechanisms'/outer/(source+'.npz');receipt=m.read(m.A/'mechanisms'/outer/(source+'.json'));assert m.sha(path)==receipt['sha256']
  with np.load(path,allow_pickle=False) as z:state=z['state'];rows=z['rows']
  lookup=np.full(int(rows.max())+1,-1,int);lookup[rows]=np.arange(len(rows))
  with np.load(m.ROOT/'experiments/v52_joint_mean_rank/analysis/labels'/(source+'.npz'),allow_pickle=False) as z:labels=z['mean']
  nmissing=0;ntotal=0;maxerr=0.
  for gi,group in enumerate(m.source_meta(source)['groups']):
   ti=group['target_index']
   for view in range(2):
    total=0;actual=np.zeros(len(genes));control=np.zeros(len(genes))
    for lane in group['lanes']:
     tr=lookup[np.asarray(lane['treated_views'][view],int)];cr=lookup[np.asarray(lane['control_views'][view],int)];nt=np.bincount(state[tr],minlength=4);nc=np.bincount(state[cr],minlength=4);nmissing+=int(nt[nc==0].sum());ntotal+=len(tr)
     if gi<2:
      # Direct sparse population summation, without calling aggregate or group_summary.
      actual+=x[tr].sum(0).A1;control+=x[cr].mean(0).A1*len(tr);total+=len(tr)
    if gi<2:
     expected=np.clip(np.log2((actual/total+.1)/(control/total+.1)),-3,3);error=float(abs(expected-labels[ti,view]).max());assert error<2e-6;maxerr=max(maxerr,error);checks.append(dict(source=source,target=group['target'],view=view,max_absolute_LFC_label_error=error))
  empty.append(dict(source=source,treated_state_cells_without_same_lane_state_controls=nmissing,total_treated_view_cells=ntotal,fraction=nmissing/max(ntotal,1)));print('V61 population labels independently replayed',source,maxerr,flush=True)
 m.write(m.A/'population_statistics_audit.json',dict(passed=True,checks=checks,empty_lane_state_accounting=empty,auditor_sha256=m.sha(__file__),plan_sha256=m.sha(m.A/'plan.json'),scope='Replays two full matched target populations in both disjoint views per source; counts empty lane-state imputation exposure across all eligible targets. Does not establish causal state identity.',D3_opened=False))

if __name__=='__main__':
 with v.m.threadpool_limits(limits=2):main()
