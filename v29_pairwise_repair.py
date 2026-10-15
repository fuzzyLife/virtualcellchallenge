"""Pair-aware integer CPM correction; exact libraries and gene totals preserved.

Each accepted 2x2 swap strictly reduces the combined absolute CPM residual of
its two genes, hence the global L1 objective. No increase in the5% movement
budget. This is a numerical solver, not a treatment-dependent calibration fit.
"""
import numpy as np
import heapq
from v29_swap_kernel import best_swap

def repair_pairwise(count,reference,max_fraction=.05,tolerance=.0001,max_steps=20000):
 x=np.asarray(count,np.int32).copy();ref=np.asarray(reference,np.int32)
 libs=ref.sum(1,dtype=np.int64);inv=1/libs;target=ref.T@inv;residual=target-x.T@inv
 initial=float(abs(residual).sum()/len(x));budget=int(libs.sum()*max_fraction);moved=0;steps=0
 stamps=np.zeros(len(residual),np.int64)
 def heaps():
  positive=[(-float(residual[g]),int(g),int(stamps[g])) for g in np.flatnonzero(residual>1e-10)]
  negative=[(float(residual[g]),int(g),int(stamps[g])) for g in np.flatnonzero(residual<-1e-10)]
  heapq.heapify(positive);heapq.heapify(negative);return positive,negative
 hp,hn=heaps()
 def top(heap,width):
  chosen=[]
  while heap and len(chosen)<width:
   item=heapq.heappop(heap)
   if item[2]==stamps[item[1]]:chosen.append(item)
  for item in chosen:heapq.heappush(heap,item)
  return [item[1] for item in chosen]
 for step in range(max_steps):
  if abs(residual).sum()/len(x)<=tolerance or moved+2>budget:break
  pos=top(hp,16);neg=top(hn,16)
  if not pos or not neg:break
  best=None
  for width in [4,16]:
   result=best_swap(x,inv,residual,np.asarray(pos[:width],np.int64),np.asarray(neg[:width],np.int64),budget-moved)
   if result[0]>=0:best=result;break
  if best is None:break
  g,a,donor,receiver,k,shift=best;x[donor,g]-=k;x[receiver,g]+=k;x[donor,a]+=k;x[receiver,a]-=k;residual[g]-=shift;residual[a]+=shift;moved+=2*k;steps+=1
  for gene in [g,a]:
   stamps[gene]+=1
   if residual[gene]>1e-10:heapq.heappush(hp,(-float(residual[gene]),gene,int(stamps[gene])))
   elif residual[gene]<-1e-10:heapq.heappush(hn,(float(residual[gene]),gene,int(stamps[gene])))
  if steps%500==0:residual=target-x.T@inv;stamps.fill(0);hp,hn=heaps()
 np.testing.assert_array_equal(x.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(x.sum(0,dtype=np.int64),ref.sum(0,dtype=np.int64));assert x.min()>=0
 final=float(abs(target-x.T@inv).sum()/len(x));assert final<=initial+1e-10
 return x,{'CPM_repair_molecules':moved,'CPM_repair_steps':steps,'initial_CPM_L1':initial,'mean_reference_L1':final,'CPM_repair_budget':budget,'CPM_target_tolerance':tolerance,'CPM_budget_exhausted':moved+2>budget}
