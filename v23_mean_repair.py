"""V23 integer population surgery: match V16 gene totals and CPM means."""
import numpy as np

def preserve_reference_means(proposed,reference,seed=23,max_steps=20000):
 assert np.issubdtype(np.asarray(proposed).dtype,np.integer) and np.issubdtype(np.asarray(reference).dtype,np.integer)
 x=np.asarray(proposed,dtype=np.int32).copy();ref=np.asarray(reference,dtype=np.int32)
 assert x.shape==ref.shape and (x>=0).all() and (ref>=0).all()
 libs=x.sum(1,dtype=np.int64);assert (libs>0).all();np.testing.assert_array_equal(libs,ref.sum(1,dtype=np.int64));inv=1/libs.astype(float);target=ref.T@inv;totals=ref.sum(0,dtype=np.int64);delta=totals-x.sum(0,dtype=np.int64)
 initial_error=float(abs(target-x.T@inv).sum()/len(x));moved_bulk=0
 negative=list(np.flatnonzero(delta<0));positive=list(np.flatnonzero(delta>0));i=j=0
 while i<len(negative) and j<len(positive):
  a=negative[i];b=positive[j];need=min(int(-delta[a]),int(delta[b]));rows=np.flatnonzero(x[:,a]>0);amount=np.minimum(x[rows,a].astype(np.int64),np.maximum(need-np.r_[0,np.cumsum(x[rows,a].astype(np.int64))[:-1]],0));x[rows,a]-=amount.astype(np.int32);x[rows,b]+=amount.astype(np.int32);delta[a]+=need;delta[b]-=need;moved_bulk+=need
  if delta[a]==0:i+=1
  if delta[b]==0:j+=1
 np.testing.assert_array_equal(x.sum(0,dtype=np.int64),totals)
 residual=target-x.T@inv;weights=1/np.maximum(target,.002);blocked=np.zeros(x.shape[1],bool);rng=np.random.default_rng(seed);molecules=0;steps=0;budget=min(200000,int(libs.sum()*.02))
 def sample(rows):
  if len(rows)<=32:return rows
  order=np.argsort(inv[rows]);return np.unique(np.r_[rows[order[:8]],rows[order[-8:]],rng.choice(rows,16,replace=False)])
 for step in range(max_steps):
  pos=np.flatnonzero((residual>2e-7)&~blocked);neg=np.flatnonzero(residual< -2e-7)
  if len(pos)==0 or len(neg)==0 or molecules>=budget:break
  g=int(pos[np.argmax(residual[pos]*weights[pos])]);partners=neg[np.argsort(residual[neg]*weights[neg])[:8]];best=None;gain_best=0.
  high=sample(np.flatnonzero(x[:,g]>0))
  for k in partners:
   low=sample(np.flatnonzero(x[:,k]>0))
   if not len(low) or not len(high):continue
   d=inv[low,None]-inv[high][None,:];wish=min(residual[g],-residual[k]);amount=np.divide(wish,d,out=np.zeros_like(d),where=d>0);amount=np.rint(amount).astype(np.int64);amount=np.minimum(amount,np.minimum(x[low,k][:,None],x[high,g][None,:]));amount=np.minimum(amount,(budget-molecules)//2);amount=np.maximum(amount,0);shift=amount*d
   gain=(abs(residual[g])-abs(residual[g]-shift))*weights[g]+(abs(residual[k])-abs(residual[k]+shift))*weights[k];gain[amount<=0]=0;idx=np.unravel_index(np.argmax(gain),gain.shape)
   if gain[idx]>gain_best:
    gain_best=float(gain[idx]);best=(int(low[idx[0]]),int(high[idx[1]]),int(k),int(amount[idx]),float(shift[idx]))
  if best is None:blocked[g]=True;continue
  a,b,k,amount,shift=best;assert a!=b and amount>0;x[a,g]+=amount;x[b,g]-=amount;x[a,k]-=amount;x[b,k]+=amount;residual[g]-=shift;residual[k]+=shift;molecules+=2*amount;steps+=1
 np.testing.assert_array_equal(x.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(x.sum(0,dtype=np.int64),totals);assert (x>=0).all();real=target-x.T@inv;np.testing.assert_allclose(real,residual,atol=1e-12)
 return x,{'mean_CPM_probability_L1_before':initial_error,'mean_CPM_probability_L1_after':float(abs(real).sum()/len(x)),'exact_V16_gene_totals':True,'exact_cell_libraries':True,'bulk_repair_molecules':moved_bulk,'CPM_cycle_molecules':molecules,'CPM_cycle_steps':steps,'relative_mean_absolute_error':float(np.mean(abs(real)*weights)),'max_mean_CPM_difference':float(abs(real).max()/len(x)*1e6)}
