"""Plot fitted source heterogeneity and paired development metric changes."""
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v17_methods import A,RAW_FAMILIES

def main():
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for f in RAW_FAMILIES:
  d=pd.read_csv(A/f'{f}_variation_fits.csv');v=d.groupby('target').cv2.mean();axes[0].hist(v,bins=np.linspace(0,.35,22),alpha=.45,label=f+' ('+str(len(v))+' targets)')
 axes[0].set(xlabel='Estimated response-strength variance (CV squared)',ylabel='Targets',title='Cross-fitted raw-source estimates');axes[0].legend(fontsize=8)
 with np.load(A/'models/all.npz',allow_pickle=False) as z:
  order=np.argsort(z['cv2'])[-20:];axes[1].barh(z['targets'][order],z['cv2'][order],color='#26828e')
 axes[1].set(xlabel='Pooled, shrunken response-strength variance',title='Twenty strongest supported target estimates');fig.savefig(A/'learned_variation.png',dpi=170);plt.close(fig)
 if (A/'development_panels.csv').exists():
  frame=pd.read_csv(A/'development_panels.csv');base=frame[frame.candidate=='v16'].set_index(['family','seed']);fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
  for ax,field in zip(axes,['fidelity','jaccard','shape_wasserstein']):
   rows=[]
   for name in ['learned_half','learned_full','generic','shuffled']:
    other=frame[frame.candidate==name].set_index(['family','seed']).loc[base.index];rows.append((other[field]-base[field]).groupby('family').mean().mean())
   ax.bar(['Half','Full','Generic','Shuffle'],rows,color=['#26828e','#3b528b','#999999','#bbbbbb']);ax.axhline(0,color='black',lw=.7);ax.set(title=field+' change',ylabel='Paired family-balanced difference from V16');ax.tick_params(axis='x',rotation=20)
  fig.savefig(A/'development_changes.png',dpi=170);plt.close(fig)
 print('V17 scientific figures saved',flush=True)
if __name__=='__main__':main()
