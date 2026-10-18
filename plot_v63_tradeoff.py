"""Plot measured error/identity trade-offs with their different cohorts explicit."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v63_identity_guard/analysis'
r=json.loads((A/'count_profile_summary.json').read_text());b=json.loads((A/'broad_identity_results.json').read_text());order=['H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2'];labels=['H1','HCT116','HEK293T','CD4 rest','CD4 stimulated']
error=[100*(r['by_context'][s]['ratio_of_summed_profile_errors']-1) for s in order];change=[]
for s in order:
 z=next(v for v in b['records'] if v['family']==s);change.append(z['metrics']['guarded_assay']['pds']-z['metrics']['v24']['pds'])
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for a,v,title,xlabel in [(ax[0],error,'Generated cells: profile error','Change versus V24 (%) — lower is better'),(ax[1],change,'Desired means: broader identity','Raw PDS change versus V24 — higher is better')]:
 a.barh(labels,v,color=['#267b68' if (x<0 if a is ax[0] else x>0) else '#b66950' for x in v]);a.axvline(0,color='#4b5563',lw=.9);a.set_title(title,loc='left',fontweight='bold');a.set_xlabel(xlabel);a.invert_yaxis();a.spines[['top','right','left']].set_visible(False);a.tick_params(axis='y',length=0);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
fig.suptitle('V63 trade-off: improved expression error does not establish an overall gain',fontsize=12,fontweight='bold')
fig.text(.5,-.035,'Different public cohorts: 30 target/context combinations for cell errors; 655 target/context pairs for mean identity. No private score.',ha='center',fontsize=9)
path=A/'tradeoff.png';fig.savefig(path,dpi=150,bbox_inches='tight',facecolor='white');plt.close(fig);print(path,flush=True)
