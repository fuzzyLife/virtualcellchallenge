"""Post-hoc correction of V59/V60 oracle aggregation; no model refit or selection."""
import argparse
import numpy as np
import pandas as pd
import v59_state_composition as m
from v61_pooled_state_response import aggregate

def main(version):
 if version==59:
  from run_v59_transfer_stage import activate
 else:
  from v60_depth_adjusted_states import activate
 activate();records=[];genes,targets=m.axes();nonpanel=~np.isin(genes,targets)
 for source in m.SOURCES:
  x,meta=m.source_data(source);outer=m.GROUP[source];file=m.CACHE/'mechanisms'/outer/(source+'.npz');assert m.sha(file)==m.read(m.A/'mechanisms'/outer/(source+'.json'))['sha256']
  with np.load(file,allow_pickle=False) as z:state=z['state'];rows=z['rows']
  np.testing.assert_array_equal(rows,meta['rows']);lookup=np.full(int(rows.max())+1,-1,int);lookup[rows]=np.arange(len(rows))
  for group in m.source_meta(source)['groups']:
   own=list(genes).index(group['target']);lanes=[]
   for lane in group['lanes']:
    vv=[]
    for view in range(2):
     tr=lookup[np.asarray(lane['treated_views'][view],int)];cr=lookup[np.asarray(lane['control_views'][view],int)];nt,tm,t0=m.group_summary(x,state,tr);nc,cm,c0=m.group_summary(x,state,cr);tilt,within,valid=m.estimate(nt,nc,tm,cm);valid[:,~meta['measured']]=False;valid[:,own]=False;within[~valid]=0;vv.append(dict(raw=(nt,nc,tm,cm,t0,c0),tilt=tilt,within=within))
    lanes.append(vv)
   pooled=[aggregate([lane[view]['raw'] for lane in lanes]) for view in range(2)]
   for fit_view in range(2):
    test=pooled[1-fit_view];train=pooled[fit_view];total=test['n'];pred={k:np.zeros(len(genes)) for k in ['composition','within','joint','state_permuted']}
    for pair in lanes:
     tr=pair[fit_view];nt,nc,tm,cm,t0,c0=pair[1-fit_view]['raw'];weight=nt.sum()/total
     pred['composition']+=weight*m.reconstruct(nc,cm,tr['tilt']);pred['within']+=weight*m.reconstruct(nc,cm,np.zeros(4),tr['within']);pred['joint']+=weight*m.reconstruct(nc,cm,tr['tilt'],tr['within']);pred['state_permuted']+=weight*m.reconstruct(nc,cm,np.roll(tr['tilt'],1))
    pred={k:m.logratio(value,test['control']) for k,value in pred.items()};pred['zero']=np.zeros(len(genes));pred['direct_mean']=m.logratio(train['actual'],train['control']);truth=m.logratio(test['actual'],test['control']);axis=meta['measured']&nonpanel&(test['control']>5);axis[own]=False
    for arm,value in pred.items():records.append(dict(source=source,target=group['target'],fit_view=fit_view,arm=arm,MSE=float(np.mean((value[axis]-truth[axis])**2))))
  print('Population diagnostic corrected',version,source,flush=True)
 frame=pd.DataFrame(records);frame.to_csv(m.A/'population_oracle_cases.csv',index=False);summary=frame.groupby(['source','arm']).MSE.mean().unstack();summary.to_csv(m.A/'population_oracle_summary.csv')
 m.write(m.A/'population_oracle_correction.json',dict(complete=True,post_hoc_diagnostic=True,reason='Original lane-level errors emphasized underpowered tiny lanes. Aggregate reconstructed CPM and matched actual/control populations before response error. Existing predictions, fitted mechanisms, thresholds and original oracle files are unchanged.',version=version,script_sha256=m.sha(__file__),plan_sha256=m.sha(m.A/'plan.json'),cases_sha256=m.sha(m.A/'population_oracle_cases.csv'),D3_opened=False))
 print(summary.to_string(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--version',type=int,choices=[59,60],required=True);a=p.parse_args()
 with m.threadpool_limits(limits=2):main(a.version)
