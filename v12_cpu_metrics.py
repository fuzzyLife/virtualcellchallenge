"""Exact Scanpy CPU normalization/log-rank arithmetic; isolated from stage-1 screen.
Stage 1 uses mathematically equivalent CPM ranks, which can differ at floating ties.
"""
import numpy as np
from scipy import sparse,stats
from v12_methods import direction_metrics

def de(x,control,genes):
 x=sparse.csr_matrix(x);control=sparse.csr_matrix(control)
 xd=np.asarray(x.sum(axis=1)).ravel();cd=np.asarray(control.sum(axis=1)).ravel()
 assert (xd>0).all() and (cd>0).all()
 xn=x.astype(np.float64);cn=control.astype(np.float64)
 xn.data=np.divide(xn.data,np.repeat(xd/1e6,np.diff(xn.indptr)))
 cn.data=np.divide(cn.data,np.repeat(cd/1e6,np.diff(cn.indptr)))
 xm=np.asarray(xn.mean(axis=0)).ravel();cm=np.asarray(cn.mean(axis=0)).ravel()
 xl=xn.copy();cl=cn.copy();xl.data=np.log1p(xl.data);cl.data=np.log1p(cl.data)
 eligible=cm>5;cols=np.flatnonzero(eligible);p=np.ones(len(genes));q=p.copy()
 for start in range(0,len(cols),256):
  c=cols[start:start+256]
  p[c]=stats.mannwhitneyu(xl[:,c].toarray(),cl[:,c].toarray(),axis=0,method='asymptotic',use_continuity=False).pvalue
 p=np.nan_to_num(p,nan=1.);q[eligible]=stats.false_discovery_control(p[eligible],method='bh')
 return {'p':p,'q':q,'lfc':np.log2((xm+1e-9)/(cm+1e-9)),'eligible':eligible}
