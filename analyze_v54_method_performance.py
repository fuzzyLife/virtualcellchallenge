"""Retrospective method/leaderboard associations; no training, upload, or API calls.

Run in the project's scientific Python environment. --check verifies archived
numbers, source hashes and documentation without rewriting any result.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import platform

os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/v54_method_performance'
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')
LABELS = {'pds': 'Discrimination', 'mse': 'Expression accuracy',
          'nmae': 'Fold-change accuracy', 'fid': 'Direction fidelity',
          'reach': 'Direction reach', 'jac': 'DE overlap'}
# Explicitly encode the RELEASED arm, not every idea tested in that directory.
# Comparison versions are scientific references, not necessarily code ancestors.
METHODS = [
 (1, None, None, 'Public K562 response transfer', 'Initial public-transfer baseline', 'build_submission.py'),
 (2, 'v2_cd4_hybrid', 1, 'Donor-calibrated CD4 hybrid', 'Bundled donor/source revision', 'PROJECT_RETROSPECTIVE.md'),
 (4, 'v4_lineage_response', 2, 'Lineage-guided response coverage', 'Bundled source/coverage revision', 'V4_MODEL_STRATEGY.md'),
 (5, 'v5_atlas_calibrated', 4, 'Atlas calibration and dual moments', 'Bundled calibration/emission change', 'V5_MODEL_STRATEGY.md'),
 (6, 'v6_fidelity_calibrated', 4, 'Fidelity-guarded lineage response', 'Calibration relative to V4 benchmark', 'V6_MODEL_STRATEGY.md'),
 (7, 'v7_multisource_adaptive', 4, 'V4 effects with fallback resampling', 'V4 response effects retained; sampling changed', 'V7_LEADERBOARD_RESULTS.md'),
 (9, 'v9_activation_preserving', 4, 'V4 effects with activation-preserving emission', 'V4 response effects retained; emission changed', 'V9_LEADERBOARD_RESULTS.md'),
 (10, 'v10_population_transport', 4, 'Pooled interventions and compact population transport', 'Bundled response/dose/shrinkage/emission revision versus previous best', 'V10_STRATEGY.md'),
 (11, 'v11_consensus_coverage', 10, 'Consensus coverage with noise 0.65', 'consensus20_noise65: response AND noise changed', 'V11_STRATEGY.md'),
 (13, 'v13_response_programs', 11, 'Assay-aware unknown-knockdown calibration', 'Selected unknown_kd arm; response programs not selected', 'V13_STRATEGY.md'),
 (16, 'v16_control_preserving', 13, 'V13 response with residual10 emission', 'Same production response; 10% empirical residual emitter', 'V16_STRATEGY.md'),
 (20, 'v20_stack_transfer', 16, 'Bounded Stack mean expert', 'stack_mean with inherited residual10 emitter', 'V20_STRATEGY.md'),
 (24, 'v24_cd4_response', 16, 'Context-matched CD4 response expert', 'cdf50_de with inherited residual10 emitter', 'V24_STRATEGY.md'),
 (36, 'v36_population_views', 24, 'Quarter-strength retrieval contrast', 'retrieval_contrast25; full population-view model rejected', 'V35_V36_STRATEGY.md'),
 (50, 'v50_response_prior', 24, 'Generic randomized response-amplitude prior', 'randomized_prior; signs and zeros retained, no target-specific claim', 'V45_V50_RESULTS.md'),
]
# These overlap and are confounded with chronology; this is descriptive coding.
FEATURES = {
 'V10 pooled-response/compact-emission family': {10,11,13,16,20,24,36,50},
 'V11 consensus-coverage family': {11,13,16,20,24,36,50},
 'V16 residual10 emission': {16,20,24,36,50},
 'V24 context-matched CD4 expert': {24,36,50},
 'Target-specific learned residual (Stack or retrieval)': {20,36},
 'V50 generic randomized amplitude prior': {50},
 'V5 atlas/dual-moment revision': {5},
}
INPUTS = {}

def sha(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()

def read(name):
 p=ROOT/name
 INPUTS[name]=sha(p)
 return json.loads(p.read_text(encoding='utf-8-sig'))

def touch(name):
 INPUTS[name]=sha(ROOT/name)

def dump(d):
 return (json.dumps(d,indent=2,allow_nan=False)+'\n').encode('utf-8')

def csv(df):
 return df.to_csv(index=False,lineterminator='\n',float_format='%.17g').encode('utf-8')

def corr(x,y):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0: return (None,None)
 return float(pearsonr(x,y).statistic),float(spearmanr(x,y).statistic)

def loo(x,y):
 vals=[corr(np.delete(x,i),np.delete(y,i))[1] for i in range(len(x))]
 vals=[v for v in vals if v is not None]
 return (min(vals),max(vals)) if vals else (None,None)

def history():
 rows=[]
 for version,directory,base,method,scope,evidence in METHODS:
  name='analysis/leaderboard_result.json' if directory is None else f'experiments/{directory}/analysis/'+('published_result.json' if version==50 else 'leaderboard_result.json')
  d=read(name); scores=d['scaled_scores'] if version==50 else {'overall':d['overall'],**d['metrics']}
  assert (d.get('final_status') or d.get('status'))=='published'
  assert abs(sum(scores[k] for k in METRICS)/6-scores['overall'])<1e-12
  touch(evidence)
  rows.append({'version':version,'comparison_version':base,'method':method,'comparison_scope':scope,
   'entry':d.get('entry_id',d.get('entry')),'model_name':d.get('model_name',d.get('model')),
   **{k:float(scores[k]) for k in ('overall',*METRICS)},
   'partition':d['partition'],'panel':d['panel'],'anchor_set':d['anchor_set'],
   'provenance':d.get('source',d.get('provenance')),'receipt':name,'method_evidence':evidence})
 df=pd.DataFrame(rows)
 assert len(df)==15 and df.entry.nunique()==15
 assert all(df[k].nunique()==1 for k in ('partition','panel','anchor_set'))
 assert df.version.is_monotonic_increasing
 df['chronological_delta']=df.overall.diff()
 df['best_to_date']=df.overall.cummax()
 df['comparison_delta']=[np.nan if pd.isna(r.comparison_version) else r.overall-float(df.loc[df.version==r.comparison_version,'overall'].iloc[0]) for r in df.itertuples()]
 for feature,versions in FEATURES.items(): df[feature]=df.version.isin(versions).astype(int)
 # Prevent directory titles from silently substituting unsubmitted arms.
 assert read('experiments/v7_multisource_adaptive/analysis/leaderboard_result.json')['response_effect_arrays_equal_v4']
 for directory,expected in [('v11_consensus_coverage','consensus20_noise65'),('v13_response_programs','unknown_kd'),('v16_control_preserving','residual10'),('v20_stack_transfer','stack_mean'),('v24_cd4_response','cdf50_de')]:
  assert read(f'experiments/{directory}/analysis/reference_model.json')['candidate']==expected
 assert read('experiments/v36_population_views/analysis/selection.json')['selected_arm']=='retrieval_contrast25'
 assert read('experiments/v50_response_prior/analysis/selection.json')['selected']=='randomized_prior'
 assert read('experiments/v50_response_prior/analysis/reference_model.json')['selected_arm']=='randomized_prior'
 return df

def contributions(df):
 rows=[]
 for row in df.itertuples():
  if pd.isna(row.comparison_version): continue
  base=df.loc[df.version==row.comparison_version].iloc[0]
  vals={k:(getattr(row,k)-base[k])/6 for k in METRICS}
  assert abs(sum(vals.values())-row.comparison_delta)<1e-12
  rows.append({'version':row.version,'comparison_version':int(row.comparison_version),
   'method':row.method,'overall_delta':row.comparison_delta,**{k+'_overall_contribution':v for k,v in vals.items()}})
 return pd.DataFrame(rows)

def correlations(df):
 methods=[];components=[]
 for era,sub in [('All published',df),('V10 onward',df[df.version>=10]),('V11 onward',df[df.version>=11])]:
  for feature in FEATURES:
   x=sub[feature].to_numpy();y=sub.overall.to_numpy()
   a,b=corr(x,y);low,high=loo(x,y)
   yes=sub.loc[sub[feature]==1,'overall'];no=sub.loc[sub[feature]==0,'overall']
   methods.append({'era':era,'feature':feature,'n':len(sub),'with_method':len(yes),'without_method':len(no),
    'mean_with':float(yes.mean()) if len(yes) else np.nan,'mean_without':float(no.mean()) if len(no) else np.nan,
    'mean_difference':float(yes.mean()-no.mean()) if len(yes) and len(no) else np.nan,
    'pearson_binary':a,'spearman':b,'leave_one_out_spearman_min':low,'leave_one_out_spearman_max':high,
    'interpretation':'Descriptive; nested methods, chronology and selection confounded; no isolated causal effect'})
  for k in METRICS:
   a,b=corr(sub[k],sub.overall);low,high=loo(sub[k].to_numpy(),sub.overall.to_numpy())
   components.append({'era':era,'metric':k,'n':len(sub),'pearson':a,'spearman':b,
    'leave_one_out_spearman_min':low,'leave_one_out_spearman_max':high,
    'interpretation':'Part-whole association: overall is the mean of these six metrics; not predictive validation'})
 return pd.DataFrame(methods),pd.DataFrame(components)

def validation(df):
 comparisons=[]
 # Matched candidate and baseline within each saved public test. Not common panels.
 for version,directory,arm,base in [(16,'v16_control_preserving','residual10',13),(20,'v20_stack_transfer','stack_mean',16),(24,'v24_cd4_response','cdf50_de',16)]:
  path=f'experiments/{directory}/analysis/confirmation_result.json';d=read(path)
  assert d['candidate']==arm and d['new_candidate_passes']
  changes=next(c['changes'] for c in d['comparisons'] if c['candidate']==arm)
  comparisons.append((version,base,path,changes,version!=24,'Confirmation; public biology reused. '+('V16 mechanism refitted without Stim48hr; not exact published V16.' if version==24 else 'Missing NMAE panels; available metrics only.')))
 d=read('experiments/v36_population_views/analysis/selection.json');a=d['count_summaries']['retrieval_contrast25'];b=d['count_summaries']['v24']
 changes={k:a[k]-b[k] for k in ('pds','fidelity','reach','jaccard','nmae')}
 changes['relative_expression_worsening']=a['common_expression_ratio']/b['common_expression_ratio']-1
 comparisons.append((36,24,'experiments/v36_population_views/analysis/selection.json',changes,True,'Exploratory reused development; strict guards failed; no new confirmation.'))
 path='experiments/v50_response_prior/analysis/confirmation_summary.json';d=read(path);a=d['randomized_prior'];b=d['v24']
 changes={k:a[k]-b[k] for k in ('pds','fidelity','reach','jaccard','nmae')}
 changes['relative_expression_worsening']=a['common_expression_ratio']/b['common_expression_ratio']-1
 comparisons.append((50,24,path,changes,True,'Additional-count confirmation; reused biology and seed-sensitive.'))
 mapping={'pds':('pds',1),'mse':('relative_expression_worsening',-1),'nmae':('nmae',-1),'fid':('fidelity',1),'reach':('reach',1),'jac':('jaccard',1)}
 rows=[]
 for ver,base,path,changes,primary,scope in comparisons:
  actual=df.loc[df.version==ver].iloc[0];old=df.loc[df.version==base].iloc[0]
  for metric,(local,sign) in mapping.items():
   ld=float(changes[local])*sign;od=float(actual[metric]-old[metric])
   rows.append({'version':ver,'comparison_version':base,'metric':metric,'local_oriented_delta':ld,
    'official_scaled_delta':od,'local_improves':ld>1e-12,'official_improves':od>1e-12,
    'nonzero_sign_agreement':None if abs(ld)<1e-12 or abs(od)<1e-12 else bool(np.sign(ld)==np.sign(od)),
    'primary_diagnostic':primary,'scope':scope,'evidence':path})
 val=pd.DataFrame(rows);cor=[]
 for scope,subset in [('Four less-mismatched comparisons (V16,V20,V36,V50)',val[val.primary_diagnostic]),('Sensitivity including refitted V24 comparator',val)]:
  for metric in METRICS:
   q=subset[subset.metric==metric];x=q.local_oriented_delta.to_numpy();y=q.official_scaled_delta.to_numpy()
   a,b=corr(x,y);lo,hi=loo(x,y);nonzero=q.nonzero_sign_agreement.dropna()
   cor.append({'scope':scope,'metric':metric,'n':len(q),'pearson':a,'spearman':b,
    'leave_one_out_spearman_min':lo,'leave_one_out_spearman_max':hi,
    'nonzero_pairs':len(nonzero),'sign_agreements':int(nonzero.sum()),
    'interpretation':'Tiny selected retrospective sample; different public panels/anchors and weighting; neither calibrated nor prospective'})
 return val,pd.DataFrame(cor)

def mdtable(df,columns):
 lines=['| '+' | '.join(columns)+' |','|'+'|'.join('---' for _ in columns)+'|']
 for _,r in df.iterrows():
  vals=[]
  for k in columns:
   v=r[k]
   vals.append(f'{v:.6f}' if isinstance(v,(float,np.floating)) and pd.notna(v) else ('unavailable' if pd.isna(v) else str(v)))
  lines.append('| '+' | '.join(vals)+' |')
 return '\n'.join(lines)

def report(df,delta,mc,cc,val,vc,sens):
 best=df.loc[df.overall.idxmax()];late=cc[cc.era=='V11 onward']
 localreach=val[(val.metric=='reach')&val.primary_diagnostic]
 lines=['# V54: which methods tracked published performance?', '',
  'Retrospective analysis of 15 published entries, all on the same validation panel and anchor set. Score receipts are the source of truth; ranks vary with other teams and are not used as the outcome. No new model, predictions, upload or fresh API check is performed.', '',
  f'Best published: V{best.version}, overall **{best.overall:.9f}**. Method associations are descriptive: multiple changes travel together, variants were selected after public testing, and no randomized component ablation was uploaded.', '',
  'The V10-family indicator has Spearman rho **0.866** with overall across all 15 uploads. Residual10 has rho **0.845** within V10 onward (8 uploads). These are inherited-pipeline associations, not isolated component effects. For example, learned residuals have a positive all-history association even though both uploaded residual revisions (V20 and V36) lose to their own baselines.', '',
  '## What actually improved', '',
  mdtable(delta[['comparison_version','version','method','overall_delta']],['comparison_version','version','method','overall_delta']), '',
  'The largest benchmark gains are V4 versus V2 (+0.172652) and V10 versus V4 (+0.167428). In the V4-to-V10 change, fidelity contributes +0.110253 and PDS +0.052860 to overall: together 97.4% of the net gain. This is an arithmetic decomposition, not proof of which code component caused it. V10 versus V4 is a bundled benchmark comparison, not an isolated component test. V11 also changes both consensus coverage and noise. V16 is a clearer emission comparison because its production response remains V13; V24 changes the mean expert while retaining the V16 emitter. V20, V36 and V50 use their actual scientific baseline, rather than the immediately preceding upload.', '',
  'V7 released V4 effects with resampling. V13 selected unknown-knockdown handling, not the tested response-program learner. V36 released retrieval_contrast25; its full population-view model was rejected. V50 selected the fixed randomized amplitude prior, not demonstrated target-specific KOLF transfer. V52/V53 joint mean/rank training has no official score and is excluded. Directory names and experiment counts are not evidence that a method improved the leaderboard.', '',
  '## Statistical associations, with their limits', '',
  'Binary method features explicitly encode inherited components. Their Spearman associations with overall are below. They overlap and are confounded with chronology. A feature present only once has no replication; a constant feature within an era has no defined correlation. Do not read these values as effect sizes or probabilities of success.', '',
  mdtable(mc[mc.era.isin(['All published','V10 onward'])],['era','feature','with_method','without_method','spearman']), '',
  'The all-history correlations mainly separate the early negative-score pipeline from the later population-emission family. Within the later family, there is much less evidence to distinguish individual components. No high-dimensional regression, causal estimate or independence-based p-value is justified by 15 selected, related uploads.', '',
  'Official metric-to-overall Spearman correlations in the V11+ era (7 uploads) are:', '',
  mdtable(late,['metric','spearman','leave_one_out_spearman_min','leave_one_out_spearman_max']), '',
  '**Part-whole coupling:** overall is exactly the mean of the six scaled components. Correlation of one component with overall is therefore partly arithmetic, not an independently validated predictor. Expression and fold-change columns here are scaled accuracies: higher is better. Raw public MSE/NMAE are errors: lower is better.', '',
  '## Did our local tests predict official gains?', '',
  'Extracted candidate-minus-comparator deltas from saved result JSONs, oriented so positive is favorable. The four primary diagnostic comparisons are V16, V20, V36 and V50; V24 is reported separately because its public baseline was refitted without Stim48hr and is not the exact published V16. Even the four primary comparisons differ in study panels, scoring aggregation, missing metrics, and whether they are confirmation or reused development. These are within-test matched changes, not a common independent benchmark.', '',
  mdtable(vc[vc.scope.str.startswith('Four')],['metric','n','spearman','leave_one_out_spearman_min','leave_one_out_spearman_max','sign_agreements','nonzero_pairs']), '',
  'Leave-one-out ranges are sensitivity diagnostics, not confidence intervals. With four selected comparisons, a coefficient is too unstable to choose a model. Raw public deltas and private anchored scaled deltas have different units; they are never averaged together. Full paired rows, including V24 sensitivity, are in local_validation_pairs.csv.', '',
  'The most revealing mismatch is REACH:', '',
  mdtable(localreach,['version','comparison_version','local_oriented_delta','official_scaled_delta']), '',
  f"V50 public REACH gain is {sens['full_public_raw_reach']['delta']:+.6f}; removing CD4_Rest_D1/MTHFSD (three reference DE genes) gives {sens['leave_case_out'][0]['leave_case_out_reach_delta']:+.6f}. This previously audited sensitivity explains why a positive aggregate did not establish a robust gain. It is not a reason to remove that case from official-style summaries.", '',
  'In these four selected comparisons, expression and fidelity changes each agree in direction in 4/4 cases; REACH and NMAE agree in only 1/4. This is a descriptive sign check across heterogeneous tests, not an accuracy estimate for future releases. The sample is too small and selected to rank validation metrics reliably.', '',
  'V20 also improved public confirmation REACH but lost official REACH. V36 improved public development REACH but lost it privately. V16 is the successful case: both public and official REACH increased. Local public gains have not been a reliable release signal.', '',
  '## Implication for the next model', '',
  'Preserve V24 as the measured reference. The strongest observed progress came from response pooling, consensus coverage and count-population construction, with a smaller gain from a context-compatible CD4 expert. Learned residuals and generic amplification have not yet beaten that reference. This does not establish that neural networks cannot work: only two target-specific learned residual variants were actually uploaded, and the new joint model was not.', '',
  'Evaluate future candidates on one frozen, common study/target panel with identical controls and emission seeds. Report per-study mean-response error, rank-response error, realized count metrics, null DE and case influence separately; retain missing values. Require gains to survive whole-study exclusion, generation seeds and leave-target-out sensitivity. Keep a simple linear and shuffled-target control. Select before opening new biological reserve data. Those tests address the observed validation mismatch; this retrospective cannot forecast a leaderboard gain.', '',
  '## Reproduce', '',
  'Use the existing Python 3.12 scientific environment (NumPy, pandas, SciPy, matplotlib). From Windows Command Prompt:', '',
  '~~~bat', 'cd /d L:\\promec\\Animesh\\Download\\Virtual',
  '"C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe" analyze_v54_method_performance.py',
  '"C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe" analyze_v54_method_performance.py --check', '~~~', '',
  'Outputs: [complete history](analysis/v54_method_performance/published_history.csv), [method coding](analysis/v54_method_performance/method_associations.csv), [paired component contributions](analysis/v54_method_performance/component_contributions.csv), [local validation pairs](analysis/v54_method_performance/local_validation_pairs.csv), [validation correlations](analysis/v54_method_performance/local_validation_correlations.csv), and [source hashes/audit](analysis/v54_method_performance/audit.json). --check recomputes numerical outputs and verifies all archived files/source hashes without rewriting them. The source receipts and previous scientific results remain unchanged.', '',
  '![Official score history](analysis/v54_method_performance/history.png)', '',
  '![Recent contributions and validation mismatch](analysis/v54_method_performance/diagnostics.png)', '']
 return '\n'.join(lines).encode('utf-8')

def plot(df,delta,val):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(2,1,figsize=(12,8),layout='constrained')
 for ax,sub,title in [(axs[0],df,'Published score history: same panel and anchors'),(axs[1],df[df.version>=10],'Later submissions: V24 remains highest')]:
  xx=np.arange(len(sub));yy=sub.overall.to_numpy()
  ax.plot(xx,yy,'o-',color='#236b8e',lw=1.8)
  ax.plot(xx,sub.best_to_date,'--',color='#7c7c7c',lw=1,label='Best to date')
  ax.set_xticks(xx,[f'V{v}' for v in sub.version]);ax.set_ylabel('Official overall score');ax.set_title(title,loc='left');ax.grid(axis='y',alpha=.2)
  for i,(v,y) in enumerate(zip(sub.version,yy)):
   if ax==axs[1] or v in (1,4,10): ax.annotate(f'{y:.4f}',(i,y),xytext=(0,9 if v not in (13,20,50) else -17),textcoords='offset points',ha='center',fontsize=9)
  ax.margins(y=.22)
 axs[0].axhline(0,color='black',lw=.7,alpha=.4);axs[0].legend(loc='lower right')
 fig.savefig(OUT/'history.png',dpi=155);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained',gridspec_kw={'width_ratios':[1.35,1]})
 recent=delta[delta.version.isin([16,20,24,36,50])]
 matrix=recent[[k+'_overall_contribution' for k in METRICS]].to_numpy()
 bound=float(np.max(np.abs(matrix)))
 im=axs[0].imshow(matrix,cmap='RdBu',vmin=-bound,vmax=bound,aspect='auto')
 axs[0].set_xticks(range(6),['PDS','Expr.','LFC','FID','REACH','JAC'],rotation=30)
 axs[0].set_yticks(range(len(recent)),[f'V{int(r.version)} minus V{int(r.comparison_version)}' for r in recent.itertuples()])
 axs[0].set_title('Contribution to official overall change\n(each scaled metric change divided by 6)',loc='left')
 for i in range(len(recent)):
  for j in range(6): axs[0].text(j,i,f'{matrix[i,j]:+.5f}',ha='center',va='center',fontsize=8,color='white' if abs(matrix[i,j])>bound*.65 else 'black')
 fig.colorbar(im,ax=axs[0],shrink=.7)
 q=val[(val.metric=='reach')&val.primary_diagnostic]
 axs[1].scatter(q.local_oriented_delta,q.official_scaled_delta,s=60,c='#236b8e')
 for r in q.itertuples(): axs[1].annotate(f'V{r.version}',(r.local_oriented_delta,r.official_scaled_delta),xytext=(6,7),textcoords='offset points')
 axs[1].axhline(0,color='gray',lw=.8);axs[1].axvline(0,color='gray',lw=.8)
 axs[1].set_xlabel('Public raw REACH change');axs[1].set_ylabel('Official scaled REACH change')
 axs[1].set_title('Public REACH gains often reversed\n4 selected comparisons; different public panels',loc='left');axs[1].margins(.22)
 fig.savefig(OUT/'diagnostics.png',dpi=155);plt.close(fig)

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
 touch(Path(__file__).name)
 df=history();delta=contributions(df);mc,cc=correlations(df);val,vc=validation(df)
 sens=read('analysis/v51_plateau/selection_sensitivity.json')
 for name in ('V52_JOINT_RESULTS.md','V53_JOINT_GUARD_RESULTS.md','V51_PLATEAU_DIAGNOSIS.md'):touch(name)
 tables={'published_history.csv':df,'component_contributions.csv':delta,'method_associations.csv':mc,
  'metric_associations.csv':cc,'local_validation_pairs.csv':val,'local_validation_correlations.csv':vc}
 outputs={k:csv(v) for k,v in tables.items()}
 outputs['V54_METHOD_PERFORMANCE.md']=report(df,delta,mc,cc,val,vc,sens)
 # Numeric checks guard against erroneous history, arithmetic, arms and pairing.
 assert len(val)==30 and int(df.loc[df.overall.idxmax(),'version'])==24
 assert set(val.loc[val.primary_diagnostic,'version'])=={16,20,36,50}
 assert (val.loc[(val.metric=='reach')&val.primary_diagnostic,'local_oriented_delta']>0).all()
 assert sum(val.loc[(val.metric=='reach')&val.primary_diagnostic,'official_scaled_delta']>0)==1
 audit={'scope':'Post-outcome descriptive associations; no causal or significance claim',
  'published_entries':len(df),'same_panel_and_anchors':True,'overall_mean6_max_error':float(np.max(abs(df[list(METRICS)].mean(axis=1)-df.overall))),
  'component_delta_max_error':float(np.max(abs(delta[[k+'_overall_contribution' for k in METRICS]].sum(axis=1)-delta.overall_delta))),
  'production_selected_arms_checked':True,'unsubmitted_joint_models_excluded':True,
  'best_version':24,'best_overall':float(df.overall.max()),'local_pairs':len(val),
  'primary_local_versions':[16,20,36,50],'raw_scaled_metric_directions_checked':True,
  'method_features_are_nested_and_nonrandomized':True,'leave_one_out_is_not_confidence_interval':True,
  'inputs_sha256':dict(sorted(INPUTS.items())),
  'numerical_outputs_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in outputs.items()}}
 if args.check:
  saved=json.loads((OUT/'audit.json').read_text())
  assert saved['analysis']==audit,'Analysis/input hashes changed'
  for name,body in outputs.items():
   path=ROOT/name if name.endswith('.md') else OUT/name
   assert path.read_bytes()==body,f'Archived output changed: {name}'
  for name,hashed in saved['figures_sha256'].items(): assert sha(OUT/name)==hashed
  print('PASS: 15 receipts, selected arms, shared anchors, mean-six arithmetic, comparison deltas, validation pairing, source/output hashes and figures.')
 else:
  OUT.mkdir(parents=True,exist_ok=True)
  for name,body in outputs.items(): (ROOT/name if name.endswith('.md') else OUT/name).write_bytes(body)
  plot(df,delta,val)
  saved={'analysis':audit,'runtime':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__},'figures_sha256':{k:sha(OUT/k) for k in ('history.png','diagnostics.png')}}
  (OUT/'audit.json').write_bytes(dump(saved))
  print(delta[['version','comparison_version','overall_delta']].to_string(index=False))
  print(vc[vc.scope.str.startswith('Four')][['metric','spearman','sign_agreements','nonzero_pairs']].to_string(index=False))
  print('Wrote V54_METHOD_PERFORMANCE.md and analysis/v54_method_performance; no predictions modified.')

if __name__=='__main__':main()
