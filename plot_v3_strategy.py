"""Plot source-validation evidence and the V3 context ensemble assumptions."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v3_context_ensemble/analysis'


def main():
    weights=json.loads((OUT/'context_mixture_weights.json').read_text())
    cd4=json.loads((OUT/'cd4_model_selection.json').read_text())
    h1=json.loads((OUT/'h1_model_selection.json').read_text())
    names=['K562','HCT116','HEK293T','H1','CD4']; contexts=['A','B','C']
    colors=['#aa3377','#ccbb44','#66ccee','#4477aa','#228833']
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    bottom=np.zeros(3)
    for name,color in zip(names,colors):
        values=np.array([weights[c]['cell_profiles'][name] for c in contexts])
        axes[0,0].bar(contexts,values,bottom=bottom,label=name,color=color)
        bottom+=values
    axes[0,0].set_title('Basal-expression fit weights\nThese are not cell fractions')
    axes[0,0].set_ylabel('Normalized NNLS coefficient');axes[0,0].set_ylim(0,1.08)
    axes[0,0].legend(fontsize=8,ncol=3,loc='upper center')
    bottom=np.zeros(3)
    for condition,color in zip(['Rest','Stim8hr','Stim48hr'],['#66ccee','#4477aa','#aa3377']):
        values=np.array([.5*weights[c]['cd4_conditions']['CD4_'+condition]+.5/3 for c in contexts])
        axes[0,1].bar(contexts,values,bottom=bottom,label=condition,color=color)
        bottom+=values
    axes[0,1].set_title('Regularized CD4 condition weights\nHalf fitted + half equal consensus')
    axes[0,1].set_ylabel('Nominal weight before missing-condition adjustment')
    axes[0,1].set_ylim(0,1.08);axes[0,1].legend(fontsize=8,ncol=3,loc='upper center')
    values=[cd4['v2_outer_cv_relative_zero_loss'],cd4['nested_selection_relative_zero_loss'],h1['nested_target_holdout_relative_zero_loss']]
    axes[1,0].bar(['CD4 V2','CD4 V3\nnested folds','H1 V3\nheld targets'],values,color=['#999999','#228833','#4477aa'])
    axes[1,0].axhline(1,color='black',linestyle='--',linewidth=1,label='Zero response')
    for i,value in enumerate(values):
        axes[1,0].text(i,value+.025,f'{value:.4f}',ha='center')
    axes[1,0].set_ylim(0,1.12);axes[1,0].set_ylabel('Error relative to each source zero-response anchor')
    axes[1,0].set_title('Source validation: lower is better\nCD4 and H1 use different source datasets')
    candidates=pd.read_csv(OUT/'control_template_candidates.csv')
    for context,color in zip(contexts,['#4477aa','#ee7733','#228833']):
        sub=candidates[candidates.context==context]
        axes[1,1].plot(sub.candidate,sub.control_profile_mse,'.',color=color,label=context,alpha=.6)
        best=sub.loc[sub.control_profile_mse.idxmin()]
        axes[1,1].plot(best.candidate,best.control_profile_mse,'*',color=color,markersize=14)
    axes[1,1].set_title('Balanced control-template selection\nStars mark the chosen set of 400 cells')
    axes[1,1].set_xlabel('Deterministic candidate index');axes[1,1].set_ylabel('Pooled log-profile MSE vs full control pool')
    axes[1,1].legend(title='Context')
    fig.suptitle('V3 context-specific response ensemble\nSource tests and expression fits do not estimate hidden leaderboard performance',fontsize=13)
    fig.savefig(OUT/'v3_strategy.png',dpi=180)
    fig.savefig(OUT/'v3_strategy.pdf')
    plt.close(fig)


if __name__=='__main__':
    main()
