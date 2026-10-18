"""Control-state diagnostic: shared coordinates, source identity and library effects."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v59_state_composition import *
from run_v59_closed_stage import activate

def main():
 activate();model=state_model('CD4');parts=[];records=[]
 for source in SOURCES:
  x,m=source_data(source);v,_=standardize(x,m,model['positions']);z=(v-model['pca_mean'])@model['basis'].T;ids=np.flatnonzero(m['label']=='Non-Targeting');rng=np.random.default_rng(SEED);show=rng.choice(ids,min(500,len(ids)),replace=False);state=((z[show,None,:]-model['centers'][None,:,:])**2).sum(2).argmin(1);parts.append(pd.DataFrame(dict(source=source,PC1=z[show,0],PC2=z[show,1],state=state)))
  for pc in range(2):records.append(dict(source=source,PC=pc+1,controls=len(ids),library_spearman=float(spearmanr(z[ids,pc],m['depth'][ids]).statistic),detected_genes_spearman=float(spearmanr(z[ids,pc],x[ids].getnnz(1)).statistic)))
  del x,v,z;gc.collect()
 df=pd.concat(parts,ignore_index=True);df.to_csv(A/'control_state_coordinates.csv',index=False);pd.DataFrame(records).to_csv(A/'state_quality.csv',index=False)
 fig,ax=plt.subplots(1,2,figsize=(11,4.5))
 for source,g in df.groupby('source'):ax[0].scatter(g.PC1,g.PC2,s=5,alpha=.35,label=source)
 for state,g in df.groupby('state'):ax[1].scatter(g.PC1,g.PC2,s=5,alpha=.35,label='State '+str(state))
 ax[0].set_title('Control populations: CD4-held-out basis');ax[1].set_title('Four statistical RNA states')
 for a in ax:a.set_xlabel('PC1');a.set_ylabel('PC2');a.legend(fontsize=7,markerscale=2)
 fig.tight_layout();fig.savefig(A/'control_states.png',dpi=160);plt.close(fig)
 write(A/'state_plot.json',dict(plot_sha256=sha(A/'control_states.png'),quality_sha256=sha(A/'state_quality.csv'),script_sha256=sha(Path(__file__)),interpretation='NTC-only coordinates; within-lane standardization intentionally removes source means. Overlap is not proof of biological state equivalence. Library/depth correlations are diagnostic, not newly optimized features.'))
 print('Control-state visualization complete',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
