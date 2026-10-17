"""Descriptive paired-view label noise; outcomes never enter fitting features."""
from pathlib import Path
import numpy as np
import v52_joint_response as m

def main():
 rows=[]
 for source in m.SOURCES:
  z=m.labels(source);ok=z['valid']
  for kind in ['mean','rank']:
   left=z[kind][:,0].astype(float);right=z[kind][:,1].astype(float);a=left[ok];b=right[ok]
   den=ok.sum(0);lm=np.divide(np.where(ok,left,0).sum(0),den,out=np.zeros(left.shape[1]),where=den>0);rm=np.divide(np.where(ok,right,0).sum(0),den,out=np.zeros(right.shape[1]),where=den>0)
   centered_ok=ok&(den[None,:]>=3);ca=(left-lm)[centered_ok];cb=(right-rm)[centered_ok]
   rows.append(dict(source=source,study=m.GROUP[source],label=kind,eligible_pairs=len(a),
    view_correlation=float(np.corrcoef(a,b)[0,1]),mean_squared_view_disagreement=float(np.mean((a-b)**2)),
    average_view_noise_proxy=float(np.mean((a-b)**2)/4),cross_view_product=float(np.mean(a*b)),
    target_centered_view_correlation=float(np.corrcoef(ca,cb)[0,1]),target_centered_pairs=len(ca)))
 m.write(m.A/'label_noise.json',dict(records=rows,plan_sha256=m.b.sha256(m.A/'plan.json'),script_sha256=m.b.sha256(Path(__file__)),
  scope='Descriptive replicate reliability; centering uses source outcomes only for this diagnostic, never as features. Reused control cells across targets and shared studies preclude interpreting these as independent confidence bounds.'))
 print('V52 paired-view reliability recorded',flush=True)
if __name__=='__main__':main()
