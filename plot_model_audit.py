"""Export the completed model audit as a scientific figure."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/model_audit'

def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    full=pd.read_csv(OUT/'full_panel_de.csv')
    source=pd.read_csv(OUT/'source_stage_and_nulls_summary.csv')
    corr=pd.read_csv(OUT/'reference_context_correlations_common_axis.csv')
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),gridspec_kw={'width_ratios':[1,1.5,1.1]})
    contexts=['A','B','C']
    means=full.groupby('context').significant_downstream_genes.mean().reindex(contexts)
    bars=axes[0].bar(contexts,means,color=['#3478a8','#2a9d8f','#d69636'])
    for bar,context in zip(bars,contexts):
        group=full[full.context==context]
        axes[0].text(bar.get_x()+bar.get_width()/2,bar.get_height()+.035,f'{bar.get_height():.2f}',ha='center')
        axes[0].text(bar.get_x()+bar.get_width()/2,.08,f'{int((group.significant_downstream_genes==0).sum())}/300\nzero calls',ha='center',color='white',fontweight='bold')
    axes[0].set_ylim(0,.75);axes[0].set_ylabel('Mean significant downstream genes / target')
    axes[0].set_title('Submitted predictions: all 900 groups',fontsize=11)
    names=['raw','reliability_only','less_threshold','original_threshold']
    x=np.arange(4);width=.37
    for i,(kind,label,color) in enumerate([('perturbed','Perturbed source groups','#3478a8'),('non_targeting_null','Untreated source groups','#b8bdc5')]):
        data=source[source.kind==kind].set_index('stage').reindex(names)
        values=data.genes_gt_point1_at_shrinkage_point3_median.to_numpy()
        bars=axes[1].bar(x+(i-.5)*width,values,width,label=label,color=color)
        for bar,value in zip(bars,values):axes[1].text(bar.get_x()+bar.get_width()/2,value*1.12,f'{value:g}',ha='center',fontsize=8)
    axes[1].set_yscale('log');axes[1].set_ylim(.7,1800)
    axes[1].set_xticks(x,['Raw effects','Reliability\nshrinkage','Weaker filter\nz=0.5','Original filter\nz=1.96'])
    axes[1].set_ylabel('Median genes with |transferred log2 effect| > 0.1')
    axes[1].set_title('Relaxed filtering mostly adds source noise',fontsize=11)
    axes[1].legend(loc='upper right',fontsize=8,frameon=False)
    table=corr.pivot(index='source',columns='context',values='log1p_cpm_correlation').reindex(['K562','HCT116','HEK293T'])[contexts]
    heat=axes[2].imshow(table.to_numpy(),cmap='viridis',vmin=.2,vmax=.6,aspect='auto')
    for (i,j),value in np.ndenumerate(table.to_numpy()):axes[2].text(j,i,f'{value:.3f}',ha='center',va='center',color='white' if value<.45 else '#17202a')
    axes[2].set_xticks(range(3),contexts);axes[2].set_yticks(range(3),table.index)
    axes[2].set_title('Basal similarity on 7,545 common genes',fontsize=11)
    fig.colorbar(heat,ax=axes[2],fraction=.06,pad=.04,label='Correlation of log1p mean CPM')
    fig.suptitle('Published VCC baseline: near-control predictions and weak source specificity',fontsize=14,y=1.01)
    fig.text(.01,-.025,'Diagnostics use provided controls; these are not hidden-response accuracy scores. Orion uses one pilot batch per cell line.',fontsize=9,color='#555555')
    fig.tight_layout()
    fig.savefig(OUT/'audit_overview.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main()
