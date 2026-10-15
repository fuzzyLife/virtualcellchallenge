"""Pair-aware integer CPM correction; exact libraries and gene totals preserved.

Each accepted 2x2 swap strictly reduces the combined absolute CPM residual of
its two genes, hence the global L1 objective. No increase in the5% movement
budget. This is a numerical solver, not a treatment-dependent calibration fit.
"""
import numpy as np

def repair_pairwise(count,reference,max_fraction=.05,tolerance=.0001,max_steps=20000):
 x=np.asarray(count,np.int32).copy();ref=np.asarray(reference,np.int32)
 libs=ref.sum(1,dtype=np.int64);inv=1/libs;target=ref.T@inv;residual=target-x.T@inv
 initial=float(abs(residual).sum()/len(x));budget=int(libs.sum()*max_fraction);moved=0;steps=0
 for step in range(max_steps):
  if abs(residual).sum()/len(x)<=tolerance or moved+2>budget:break
  pos=np.flatnonzero(residual>1e-10);neg=np.flatnonzero(residual<-1e-10)
  if not len(pos) or not len(neg):break
  pos=pos[np.argsort(-residual[pos],kind='stable')];neg=neg[np.argsort(residual[neg],kind='stable')]
  best=None;bestgain=0.;bestcost=None
  for width in [4,16]:
   for g in pos[:width]:
    sources=np.flatnonzero(x[:,g]>0)
    if not len(sources):continue
    order=np.argsort(inv[sources],kind='stable');sources=np.unique(np.r_[sources[order[:2]],sources[np.argmax(x[sources,g])]])
    for a in neg[:width]:
     destinations=np.flatnonzero(x[:,a]>0)
     if not len(destinations):continue
     wanted=min(float(residual[g]),float(-residual[a]))
     for donor in sources:
      dst=destinations[inv[destinations]>inv[donor]]
      if not len(dst):continue
      span=float((inv[dst]-inv[donor]).max());cap=int(x[donor,g]);amounts=np.unique([1,min(cap,max(1,int(round(wanted/span)))),cap])
      for amount in amounts:
       rows=dst[x[dst,a]>=amount]
       if not len(rows):continue
       receiver=int(rows[np.argmin(abs(inv[rows]-(inv[donor]+wanted/amount)))]);d=float(inv[receiver]-inv[donor]);k=min(cap,int(x[receiver,a]),max(1,int(round(wanted/d))),(budget-moved)//2)
       if k<=0:continue
       shift=k*d;gain=abs(residual[g])+abs(residual[a])-abs(residual[g]-shift)-abs(residual[a]+shift)
       if gain>bestgain+1e-15 or (gain>1e-15 and abs(gain-bestgain)<=1e-15 and (bestcost is None or k<bestcost)):
        bestgain=float(gain);bestcost=k;best=(int(g),int(a),int(donor),receiver,k,shift)
   if best is not None:break
  if best is None:break
  g,a,donor,receiver,k,shift=best;x[donor,g]-=k;x[receiver,g]+=k;x[donor,a]+=k;x[receiver,a]-=k;residual[g]-=shift;residual[a]+=shift;moved+=2*k;steps+=1
  if steps%500==0:residual=target-x.T@inv
 np.testing.assert_array_equal(x.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(x.sum(0,dtype=np.int64),ref.sum(0,dtype=np.int64));assert x.min()>=0
 final=float(abs(target-x.T@inv).sum()/len(x));assert final<=initial+1e-10
 return x,{'CPM_repair_molecules':moved,'CPM_repair_steps':steps,'initial_CPM_L1':initial,'mean_reference_L1':final,'CPM_repair_budget':budget,'CPM_target_tolerance':tolerance,'CPM_budget_exhausted':moved+2>budget}
