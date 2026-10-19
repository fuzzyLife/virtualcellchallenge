"""Plot measured V65/V66 count trade-offs without inventing official score scales."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v65_closed_response/analysis'

def main():
 inputs={v:ROOT/('experiments/v65_closed_response/analysis/count_assessment.json' if v==65 else 'experiments/v66_assay_strength/analysis/count_assessment.json') for v in (65,66)}
 data={v:read(p) for v,p in inputs.items()};fig,axes=plt.subplots(1,2,figsize=(12.8,4.8));colors={65:'#71849c',66:'#008b7a'};studies=('H1','Orion','CD4');pos=np.arange(3)
 for i,v in enumerate((65,66)):
  r=data[v];arm=r['selected'];ratios=[r['studies'][arm][s]['profile_square_error']/r['studies']['v63'][s]['profile_square_error'] for s in studies];axes[0].vlines(pos+(i-.5)*.24,1,ratios,color=colors[v],linewidth=2);axes[0].plot(pos+(i-.5)*.24,ratios,linestyle='None',marker='o',color=colors[v],label=f'V{v} {arm}');b=r['summary']['v63'];p=r['summary'][arm]
  metrics=('pds','fidelity','reach','jaccard','nmae');gains=[100*(p[k]-b[k])/abs(b[k])*(-1 if k=='nmae' else 1) for k in metrics];gains.append(100*(1-r['paired_equal_study_profile_error_ratio']));axes[1].barh(np.arange(6)+(i-.5)*.32,gains,height=.32,color=colors[v])
 axes[0].axhline(1,color='#444444',linewidth=1,linestyle='--');axes[0].set_ylim(.978,1.008);axes[0].set_xticks(pos,studies);axes[0].set_ylabel('Generated-profile error / V63 (lower is better)');axes[0].set_title('Profile gains vary across public studies');axes[0].legend(frameon=False,fontsize=8)
 axes[1].set_yticks(np.arange(6),['PDS','Fidelity','Reach','Overlap','NMAE error','Paired profile error']);axes[1].invert_yaxis();axes[1].axvline(0,color='#444444',linewidth=1);axes[1].set_xlabel('Oriented relative raw-metric change (%)\npositive is favorable; metrics have different scales');axes[1].set_title('V66: four small gains; reach guard fails')
 fig.suptitle('V65/V66 versus V63: reused public count diagnostics, not leaderboard forecasts',fontsize=12);fig.tight_layout(rect=[0,0,.995,.94]);p=A/'count_tradeoff.png';fig.savefig(p,dpi=160);plt.close(fig)
 write(A/'count_tradeoff_receipt.json',dict(inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs.values()},figure_sha256=sha(p),plotter_sha256=sha(__file__),scope='Per-study paired profile ratios and within-metric relative raw changes. Not official anchor-scaled metric contributions; do not average the percent bars as an overall score.'))
 print('Saved the measured count trade-off plot.')
if __name__=='__main__':main()
