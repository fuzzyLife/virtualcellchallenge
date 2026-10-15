"""Exact scalar swap search compiled without fast-math or parallel reductions."""
import numpy as np
from numba import njit

@njit(cache=True,fastmath=False,nogil=True)
def best_swap(x,inv,residual,positive,negative,remaining):
 bg=-1;ba=-1;bd=-1;br=-1;bk=0;bs=0.;bestgain=0.;bestcost=9223372036854775807
 for gi in range(len(positive)):
  g=positive[gi];first=-1;second=-1;maximum=-1;maxcount=-1
  for row in range(len(x)):
   count=x[row,g]
   if count<=0:continue
   if count>maxcount:maximum=row;maxcount=count
   if first<0 or inv[row]<inv[first]:second=first;first=row
   elif second<0 or inv[row]<inv[second]:second=row
  if first<0:continue
  donors=np.sort(np.array([first,second,maximum],np.int64))
  for ai in range(len(negative)):
   a=negative[ai];wanted=min(residual[g],-residual[a])
   for di in range(3):
    donor=donors[di]
    if donor<0 or (di>0 and donor==donors[di-1]):continue
    span=0.
    for row in range(len(x)):
     if x[row,a]>0 and inv[row]>inv[donor]:span=max(span,inv[row]-inv[donor])
    if span<=0:continue
    cap=int(x[donor,g]);middle=min(cap,max(1,int(np.rint(wanted/span))));amounts=np.sort(np.array([1,middle,cap],np.int64))
    for mi in range(3):
     amount=amounts[mi]
     if mi>0 and amount==amounts[mi-1]:continue
     desired=inv[donor]+wanted/amount;receiver=-1;distance=np.inf
     for row in range(len(x)):
      if x[row,a]>=amount and inv[row]>inv[donor]:
       value=abs(inv[row]-desired)
       if value<distance:distance=value;receiver=row
     if receiver<0:continue
     d=inv[receiver]-inv[donor];k=min(cap,int(x[receiver,a]),max(1,int(np.rint(wanted/d))),remaining//2)
     if k<=0:continue
     shift=k*d;gain=abs(residual[g])+abs(residual[a])-abs(residual[g]-shift)-abs(residual[a]+shift)
     if gain>bestgain+1e-15 or (gain>1e-15 and abs(gain-bestgain)<=1e-15 and k<bestcost):
      bestgain=gain;bestcost=k;bg=g;ba=a;bd=donor;br=receiver;bk=k;bs=shift
 return bg,ba,bd,br,bk,bs
