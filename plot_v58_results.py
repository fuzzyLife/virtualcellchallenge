"""Plot V58 response errors; descriptive public-cohort comparisons, no private score."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v58_pie import *

def main():
 f=pd.read_csv(A/'response_cases.csv');f=f[f.transfer].groupby(['source','seen']).mean(numeric_only=True)
 labels=[('Rest / unseen' if s=='CD4_Rest_D2' and not seen else 'Rest / seen' if s=='CD4_Rest_D2' else 'Stim / unseen' if not seen else 'Stim / seen') for s,seen in f.index]
 x=np.arange(len(f));fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
 for ax,metric,baseline,title in [(axes[0],'mean_MSE','baseline_mean_MSE','Mean error relative to V24'),(axes[1],'rank_MSE','zero_rank_MSE','Rank error relative to zero shift')]:
  for j,(name,label,color) in enumerate([('joint','PIE heads + evidence','#267caa'),('shuffled','Trained shuffled control','#ae8c64'),('evidence','Evidence only','#61956e')]):
   y=100*(f[name+'_'+metric]/f[baseline]-1);ax.bar(x+(j-1)*.24,y,.23,label=label,color=color)
  ax.axhline(0,color='#333333',linewidth=.8);ax.set_title(title,loc='left');ax.set_ylabel('Error change (%)');ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
 fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,.925));axes[1].set_xticks(x,labels)
 fig.suptitle('CD4 transfer: lower error is better',fontsize=14);fig.text(.5,.015,'Public development cohorts; point estimates, not private leaderboard scores.',ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.04,1,.87]);path=A/'response_errors.png';fig.savefig(path,dpi=160);plt.close(fig)
 write(A/'response_plot.json',{'sha256':sha(path),'script_sha256':sha(Path(__file__)),'source_sha256':sha(A/'response_cases.csv')})
 print('Response comparison plot saved',flush=True)

if __name__=='__main__':main()
