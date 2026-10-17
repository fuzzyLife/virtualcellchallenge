"""KL projection and bounded matrix rounding; exact integer row/column margins."""
import numpy as np
from numba import njit
@njit(cache=True)
def repair_binary(up,frac,delta):
 n,g=up.shape;neg=np.where(delta<0)[0];pos=np.where(delta>0)[0];ii=0;jj=0;moves=0
 while ii<len(neg) and jj<len(pos):
  a=neg[ii];b=pos[jj];best=-1;score=-1e300
  for r in range(n):
   if up[r,a] and not up[r,b] and frac[r,b]>0:
    value=frac[r,b]-frac[r,a]
    if value>score:score=value;best=r
  if best>=0:up[best,a]=False;up[best,b]=True
  else:
   # Alternating path through at most n rows; transfers stay binary.
   prev=np.full(n,-2,np.int64);via=np.full(n,-1,np.int64);queue=np.empty(n,np.int64);head=0;tail=0;finish=-1
   for r in range(n):
    if up[r,a]:prev[r]=-1;queue[tail]=r;tail+=1
   while head<tail and finish<0:
    r=queue[head];head+=1
    if not up[r,b] and frac[r,b]>0:finish=r;break
    for rr in range(n):
     if prev[rr]!=-2:continue
     choice=-1;value=-1e300
     for c in range(g):
      if not up[r,c] and frac[r,c]>0 and up[rr,c]:
       v=frac[r,c]-frac[rr,c]
       if v>value:value=v;choice=c
     if choice>=0:prev[rr]=r;via[rr]=choice;queue[tail]=rr;tail+=1
   if finish<0:return moves,False
   r=finish;c=b
   while prev[r]>=0:
    up[r,c]=True;old=via[r];up[r,old]=False;c=old;r=prev[r]
   up[r,c]=True;up[r,a]=False
  delta[a]+=1;delta[b]-=1;moves+=1
  if delta[a]==0:ii+=1
  if delta[b]==0:jj+=1
 return moves,True
def fit_expected(initial,libs,totals,steps=500):
 e=np.asarray(initial,float).copy();libs=np.asarray(libs,np.int64);totals=np.asarray(totals,np.int64);active=totals>0
 assert libs.sum()==totals.sum() and (libs>0).all()
 e[:,~active]=0;e[:,active]=np.maximum(e[:,active],1e-10)
 for step in range(steps):
  e*=np.divide(totals,e.sum(0),out=np.zeros(len(totals)),where=e.sum(0)>0)[None,:]
  e*=libs[:,None]/e.sum(1,keepdims=True)
  error=float(abs(e.sum(0)-totals).max())
  if error<1e-6:break
 assert error<1e-4,('KL projection failed',error)
 return e,step+1
def integerize(expected,libs,totals,seed):
 e=np.asarray(expected,float);floor=np.floor(e).astype(np.int32);frac=e-floor;n,g=e.shape
 rowneed=np.asarray(libs,np.int64)-floor.sum(1,dtype=np.int64);colneed=np.asarray(totals,np.int64)-floor.sum(0,dtype=np.int64)
 assert (rowneed>=0).all() and (colneed>=0).all() and rowneed.sum()==colneed.sum()
 up=np.zeros(e.shape,np.bool_);rng=np.random.default_rng(seed)
 for r in range(n):
  k=int(rowneed[r])
  if k:
   score=frac[r]+rng.uniform(0,1e-12,g);score[frac[r]<=0]=-1
   take=np.argpartition(score,-k)[-k:];assert (frac[r,take]>0).all();up[r,take]=True
 delta=colneed-up.sum(0,dtype=np.int64);moves,ok=repair_binary(up,frac,delta);assert ok and not np.any(delta),'Binary rounding infeasible'
 out=floor+up.astype(np.int32);np.testing.assert_array_equal(out.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(out.sum(0,dtype=np.int64),totals)
 assert np.max(abs(out-e))<=1+1e-9
 return out,int(moves)

