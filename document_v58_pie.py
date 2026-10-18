"""Document completed V58 outcomes, including failures and seen/unseen strata."""
import numpy as np
import pandas as pd
from v58_pie import *

def table(frame,cols,formats=None):
 formats=formats or {};out=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
 for _,r in frame.iterrows():
  cells=[]
  for k in cols:
   v=r[k];cells.append('missing' if pd.isna(v) else (formats[k].format(v) if k in formats else str(v)))
  out.append('| '+' | '.join(cells)+' |')
 return '\n'.join(out)

def main():
 assert read(A/'final_audit.json')['passed']
 response=pd.read_csv(A/'response_cases.csv');counts=pd.read_csv(A/'count_cases.csv');paired=pd.read_csv(A/'count_paired_changes.csv')
 rr=response[response.transfer].groupby(['source','seen']).mean(numeric_only=True).reset_index()
 rr['Mean MSE gain %']=100*(1-rr.joint_mean_MSE/rr.baseline_mean_MSE)
 rr['Mean gain vs shuffled %']=100*(1-rr.joint_mean_MSE/rr.shuffled_mean_MSE)
 rr['Mean gain vs evidence-only %']=100*(1-rr.joint_mean_MSE/rr.evidence_mean_MSE)
 rr['Rank MSE gain vs zero %']=100*(1-rr.joint_rank_MSE/rr.zero_rank_MSE)
 rr['Raw PIE delta / zero MSE']=rr.raw_PIE_delta_MSE/rr.zero_delta_MSE
 rr['Calibrated delta / zero MSE']=rr.joint_delta_MSE/rr.zero_delta_MSE
 cc=counts.groupby(['source','arm']).mean(numeric_only=True).reset_index()
 pc=paired[paired.source.str.startswith('CD4')].groupby(['source','target','seen','arm']).mean(numeric_only=True).reset_index()
 paired_summary=pc.groupby(['seen','arm']).mean(numeric_only=True).reset_index()
 paired_summary.to_csv(A/'held_CD4_paired_summary.csv',index=False)
 tests={}
 for source in ['CD4_Rest_D2','CD4_Stim48_D2']:
  r=response[response.source==source].mean(numeric_only=True);c=cc[cc.source==source].set_index('arm');j=c.loc['pie_joint'];b=c.loc['v24'];s=c.loc['shuffled_joint']
  tests[source]={'mean_MSE_beats_V24':bool(r.joint_mean_MSE<r.baseline_mean_MSE),'rank_MSE_beats_zero':bool(r.joint_rank_MSE<r.zero_rank_MSE),'count_profile_MSE_beats_V24':bool(j.profile_MSE<b.profile_MSE),'FID_not_worse':bool(j.fidelity>=b.fidelity),'Jaccard_not_worse':bool(j.jaccard>=b.jaccard),'FID_beats_shuffled':bool(j.fidelity>s.fidelity),'profile_MSE_beats_shuffled':bool(j.profile_MSE<s.profile_MSE)}
 eligible=all(all(t.values()) for t in tests.values())
 leave=[]
 joint=pc[pc.arm=='pie_joint']
 for target in sorted(joint.target.unique()):
  z=joint[joint.target!=target];v=z.groupby('source').mean(numeric_only=True).mean(numeric_only=True)
  leave.append({'excluded_target':target,**{k:float(v[k]) for k in ['fidelity','reach','jaccard','nmae','profile_MSE']}})
 pd.DataFrame(leave).to_csv(A/'held_CD4_leave_target_out.csv',index=False)
 assessment={'release_supported_by_conservative_check':eligible,'checks':tests,'check_scope':'Conservative interpretation of this public development experiment, not the originally frozen selection formula or an official/private scoring guarantee. No raw local overall score is constructed.','V24_remains_published_champion':True,'plan_sha256':sha(A/'plan.json'),'D3_opened':False,'novel_independent_biological_study_claimed':False,'count_case_comparisons':len(counts),'held_CD4_count_context_target_pairs':len(joint),'response_cases':len(response)}
 write(A/'assessment.json',assessment)
 parts=['# V58 cross-dataset PIE: measured results',
 'Completed experiment: genuine saved response evidence, three retained model outputs, source-only calibration and integer-count tests. These are public development results, not a new leaderboard score. V24 remains the best published project submission (0.142997855).',
 '## What ran',
 'The hash-verified checkpoint and 17.12 GiB of training summaries were downloaded. All 64,260 saved training rows supplied the global evidence normalization. The projected evidence contains 186 records for 154 seen panel targets; 146 targets have no matching response record. Independent arithmetic passed 36 comparisons; 96 reused source/code assets (1,495,668,215 bytes) matched their previous checksums. Strict weights, finite outputs and repeated-forward checks passed.',
 'H1/Orion diagnostics exclude matching query-response records from evidence, but their studies were used to pretrain the checkpoint. They are not unseen-study tests. CD4 treatment labels do not fit or select decoder parameters. The held-study feature guard additionally passed after poisoning every CD4 response-bank row. The earlier V57 calibration path could use CD4 labels indirectly in source baseline features and is superseded. All 366 eligible response cases and '+str(len(counts))+' count comparisons were evaluated. Public biology has been used in prior project development; seeds are not independent biological replicates. D3 remains unopened.',
 '## CD4 response accuracy before count generation',
 'Positive mean/rank gains below are improvements. The raw and calibrated delta columns are error ratios: below 1 beats a zero-shift predictor. Means use equal context-target case weights within each stratum.',
 table(rr,['source','seen','Mean MSE gain %','Mean gain vs shuffled %','Mean gain vs evidence-only %','Rank MSE gain vs zero %','Raw PIE delta / zero MSE','Calibrated delta / zero MSE'],{k:'{:.3f}' for k in ['Mean MSE gain %','Mean gain vs shuffled %','Mean gain vs evidence-only %','Rank MSE gain vs zero %','Raw PIE delta / zero MSE','Calibrated delta / zero MSE']}),
 'The mean predictor is a bounded correction to the V24 public hold-study baseline. CD4 comparator arrays are bitwise unchanged; source diagnostic baselines additionally exclude CD4 under the corrected guard. Rank is tie-aware U probability minus one half; zero rank shift is a necessary simple comparator. Source-only whole-target validation selected the ridge regularization. The trained shuffled-head model and evidence-only ablation are reported in response_cases.csv and response_summary.csv, not omitted when stronger.',
 'Eligible-only raw DE-probability and conditional-LFC checks were completed before count generation; see raw_head_summary.csv. The auxiliary AP field in count_cases.csv includes all supported genes and is not used for confidence validation, because untested genes must not be treated as negatives.',
 '## Raw DE and conditional fold-change heads',
 table(pd.read_csv(A/'raw_head_summary.csv').query("source.str.startswith('CD4')",engine='python'),['source','seen','DE_prevalence','raw_pDE_AP_eligible_only','raw_PIE_LFC_DE_direction','calibrated_intended_DE_direction'],{k:'{:.4f}' for k in ['DE_prevalence','raw_pDE_AP_eligible_only','raw_PIE_LFC_DE_direction','calibrated_intended_DE_direction']}),
 'AP is evaluated only on tested/eligible genes; DE prevalence is its chance reference. Direction agreement is on actual DE genes. These are raw-head diagnostic proportions, not the official scaled fidelity score.',
 '## Integer-count results',
 'Fidelity/reach/Jaccard are raw local higher-is-better metrics. NMAE and profile MSE are errors, lower is better. Missing NMAE is retained. Actual treatment populations contain up to 400 distinct available cells, while generated populations contain 400. Each comparison uses independent NTC reference cells and the same sampling seed across arms.',
 table(cc,['source','arm','fidelity','reach','jaccard','nmae','profile_MSE'],{k:'{:.6f}' for k in ['fidelity','reach','jaccard','nmae','profile_MSE']}),
 'Seen/unseen paired changes, common-control panel diagnostics, predicted/reference DE counts, calibration-vs-transfer distinctions and leave-target-out sensitivity are retained as CSV files. The panel calculation uses one common context control and keeps expression numerators and denominators separate; no unstable ratio is averaged into a private score forecast.',
 '## Zero-intervention control checks',
 table(pd.DataFrame([{'source':read(f)['source'],**{r['arm']+'_DE_calls':r['DE_calls'] for r in read(f)['records']}} for f in sorted((A/'counts').glob('*_null.json'))]),['source','v24_DE_calls','pie_joint_DE_calls']),
 'These checks compare synthetic zero-response populations with independent real controls. Zero called DE genes here does not establish accurate perturbation responses. The joint decoder repairs a baseline artifact while losing substantial treatment direction yield.',
 '![Response errors](experiments/v58_pie_holdstudy_guard/analysis/response_errors.png)',
 '## Interpretation and limitations',
 ('The joint decoder passes the conservative public support check. This supports a further submission-development step, not a promised private gain.' if eligible else 'The joint decoder fails the conservative public support check. This experiment does not justify replacing V24 or labeling an unvalidated file a better submission.'),
 'Detailed check outcomes:\n\n```json\n'+json.dumps(tests,indent=2)+'\n```',
 'PIE p_de is a tempered, class-weighted probability output, not BH FDR. Its conditional LFC head is not supervised on every non-DE coordinate. Its mean-log-expression shift is not interchangeable with a pooled-count profile through exponentiation. The joint decoder uses an explicitly approximate Gaussian mean/rank relation, then enforces integer libraries and global moments; actual count scoring tests this approximation.',
 'One pre-existing subset limitation remains explicit: TMEM104 resolves through its SLC38A12 alias in ESM2, but that renamed entry is missing from the reused NCBI token subset. NCBI coverage is 299/300, ESM2 300/300 and perturbation-text coverage 155/300. Missing inputs are masked. No configurations were mixed or quietly changed after freezing.',
 'Further work should improve cell-context conditioning and supervised DE/direction decisions. Their contribution would need a separate frozen test. These results do not establish that the published PIE architecture fails generally; they evaluate this anonymous-context configuration and this decoder.',
 'The independent auditor runs through an immutable-JSON read cache; the frozen auditor/calculation and protocol bytes are unchanged. Two separate GPU context workers only alter scheduling, not precision, weights, samples or seeds.',
 '## Reproduce and locate artifacts',
 '[Frozen protocol](V58_PIE_CROSS_STRATEGY.md) lists commands. Start with run_v58_response.py for resumable response-stage inference/calibration once inputs/evidence are prepared; the count stage remains separate. Scientific methods are sealed by plan SHA256 `'+sha(A/'plan.json')+'`. Root scripts and compact receipts/tables are committed; large arrays are regenerable caches or exported project artifacts. No competition API upload occurs.',
 'Key files: experiments/v58_pie_holdstudy_guard/analysis/response_summary.csv, count_summary.csv, held_CD4_paired_summary.csv, held_CD4_leave_target_out.csv, panel_metrics.csv, assessment.json and final_audit.json.']
 text='\n\n'.join(parts)+'\n';path=ROOT/'V58_PIE_CROSS_RESULTS.md';path.write_text(text,encoding='utf-8')
 print(json.dumps(assessment,indent=2),flush=True)

if __name__=='__main__':main()
