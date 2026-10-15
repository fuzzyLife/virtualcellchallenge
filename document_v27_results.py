"""Document measured V27 results without changing the frozen experiment."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v27_methods import ROOT,A,ARMS,FINALIST,load_plan,write
import build_submission as b

def table(frame,cols):
 def fmt(x):
  if isinstance(x,(float,np.floating)):return f'{x:.6f}' if np.isfinite(x) else 'missing'
  return str(x)
 return '|'+ '|'.join(cols)+'|\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('|'+ '|'.join(fmt(v) for v in row)+'|' for row in frame[cols].itertuples(index=False,name=None))

def main():
 plan=load_plan();decision=json.loads((A/'release_decision.json').read_text()) if (A/'release_decision.json').exists() else None;parts=['# V27 results: control-anchored relative responses','',('**No V27 submission is released.** The fixed finalist failed prospective checks. V24 remains the best recorded published submission, overall 0.14299785521745675.' if decision else 'Execution is incomplete or a fresh-reserve/release audit is still required. No upload-ready artifact is claimed.'),'',f'Frozen plan SHA256: `{b.sha256(A/"plan.json")}`. Historical V24--V26 methods, packages and numeric evidence remain unchanged.','']
 sampling=pd.read_csv(A/'sampling_scores.csv');parts+=['## Sampling-stage ablation','',table(sampling,['stage','null_DE_calls','detection_MAE','log_variance_MAE','mean_reference_L1']),'','This spent H1 control case uses the same real reference throughout. Stages are sequential ablations, not independent studies. Intermediate stages need not meet V24 mean constraints; their drift is explicit. The final repaired empirical population exactly reproduces the frozen V26 count digest. The raw control matrix already contains measurement noise. Added sampling and count repairs can each alter its distribution; neither is a neutral empirical reference.','']
 diagnostic=A/'sampling_artifact_diagnostic.csv'
 if diagnostic.exists():
  dd=pd.read_csv(diagnostic);parts+=['The descriptive artifact check also records signed detection errors and support changes:', '', table(dd,['stage','signed_detection_error','changed_molecule_fraction_vs_raw','panel_DE_calls','other_DE_calls']), '', 'The second multinomial draw removes409,390 nonzero cell-gene entries in this case. Depth-stratum repair subsequently loses324,340 and adds27,220 entries: a net loss297,120. These are measured sequential support changes; the later V24-reference and raw-anchor rows are separate comparators, not causal continuation stages.', '']
 completed=[]
 for key in ['H1','HCT116','HEK293T','CD4']:
  n=A/f'null_{key}_scores.csv';t=A/f'treatment_{key}_panels.csv'
  if not n.exists():continue
  completed.append(key);nf=pd.read_csv(n);rows=[]
  for arm in ARMS:
   f=nf[nf.candidate==arm];rows.append({'arm':arm,'mean_calls':float(f.null_DE_calls.mean()),'max_calls':int(f.null_DE_calls.max()),'detection_MAE':float(f.detection_MAE.mean()),'log_variance_MAE':float(f.log_variance_MAE.mean()),'max_mean_L1':float(f.mean_reference_L1.max())})
  parts+=['## '+key+' control calibration','',table(pd.DataFrame(rows),['arm','mean_calls','max_calls','detection_MAE','log_variance_MAE','max_mean_L1']),'','Three control splits are cell-disjoint within each case but reuse biology across seeds. NTC fitting cells are separate from templates, actual and reference populations. Guide/batch overlap is archived; this is not guide- or batch-disjoint validation. Null relative descriptors are zero, so relative models share their corresponding anchor populations. Passing this check alone would not establish target-specific treatment accuracy.','']
  if t.exists():
   tf=pd.read_csv(t);summ=tf.groupby('candidate',sort=False)[['pds','fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']].mean().reindex(ARMS).reset_index();parts+=['## '+key+' held-family treatment diagnostic','',table(summ,['candidate','pds','fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']),'','These are raw public-development metrics, not anchor-scaled leaderboard estimates. Treatment means remain fixed to the family-excluded V24 algorithm, not the published ABC populations. All query treatment cells are excluded from fitting and inference. In the H1 pilot, all nine arms are evaluated even if null calibration fails; this diagnostic exception was prespecified before outcomes.','']
  for stage in ['null','treatment']:
   rr=A/f'{stage}_{key}_result.json'
   if rr.exists():
    result=json.loads(rr.read_text());parts+=[f'{stage.title()} checks passed: **{result["passed"]}**. Failed checks: '+(', '.join(k for k,v in result['checks'].items() if not v) or 'none')+'.','']
 models=[]
 for f in sorted((A/'models').glob('*.json')) if (A/'models').exists() else []:
  rec=json.loads(f.read_text());models.append({'model':f.stem,'epochs':len(rec['history']),'best_validation_ELBO':min(v['validation_elbo'] for v in rec['history']),'source_contexts':', '.join(rec['training_contexts'])})
 if models:parts+=['## Fitting and audit','',table(pd.DataFrame(models),['model','epochs','best_validation_ELBO','source_contexts']),'','The panel network permits measured panel outputs while masking the directly perturbed output per example. Stable state, posterior features and response basis remain nonpanel. Same source sampling, architecture, initialization and held-cell ELBO selection as V26; better likelihood is not a transfer score.','']
 audit=A/'independent_audit.json'
 if audit.exists():
  v=json.loads(audit.read_text());parts+=[f'Independent audit passed: **{v["passed"]}**. Checked {v["original_source_contexts_checked"]} source contexts, {v["original_full_axis_rows_checked"]} complete original rows, {v["fit_only_group_means_recomputed"]} fit-only group means and {v["paired_cases_checked"]} paired cases. {len(v["full_axis_replays"])} full-axis final-emitter replays matched count digests and scorer call counts. H1 source arrays match the previously audited immutable V26 training arrays.','']
 if decision:parts+=['## Release decision','',decision['reason'],'','No new package, upload or reserve opening occurred. Broader-family fitting was not executed after the H1 rejection. No diagnostic arm was promoted after the fixed finalist failed. The old champion stays available; uploading it again would duplicate an existing result.','']
 if (A/'projection_diagnostic.csv').exists():
  dd=pd.read_csv(A/'projection_diagnostic.csv');parts+=['## Post-freeze projection diagnostic', '', table(dd,['stage','DE_calls','detection_MAE','log_variance_MAE','mean_reference_L1']), '', 'This diagnostic uses the already-spent first H1 null. Removing depth strata and retaining global CPM swaps meets all first-case envelopes, but is not an eligible V27 arm and does not revise its rejection. It motivates the separately frozen V28 experiment.', '']
 if (A/'process_execution.json').exists():
  parts+=['## Operational parallelism', '', 'After null comparisons completed, only the verified original orchestration process was stopped. Completed null tables/caches/decisions were preserved. Root run_v27_with_processes.py resumes the identical frozen pipeline using four independent joblib processes and two BLAS threads each. test_v27_processes.py verifies serial/process Atlas arithmetic; the full-axis null/treatment replay audit verifies final scientific counts. analysis/orchestration_handoff.json and process_execution.json record the change. No model, threshold, seed, score or finalist changed.', '']
 parts+=['## Limits and reproduction','','One held-out source context cannot establish broad generalization. Public development biology has been reused in prior experiments; three NTC seeds do not add independent lineages. Quantile calibration may alter joint covariance and depth dependence. Source states are observational proxies, and relative Gaussian latent draws are not identified counterfactual cells. Multiplicative relative logits cannot create expression from a zero template entry before quantile mapping or moment repair. The external V25 response basis used allowed source summaries including source held-cell validation observations; validation ELBO is not wholly independent of that preprocessing, while complete H1 exclusion remains enforced. Exact mean constraints can themselves create artifacts. The learning/output-mask and calibration ideas must be judged by matched treatment comparisons, not just control reconstruction.','','Methods, inputs, release criteria and commands: [V27_STRATEGY.md](V27_STRATEGY.md). Re-run the pinned scientific runtime; frozen outcomes refuse overwrite and verified caches can resume. Run `audit_v27_experiment.py` then `document_v27_results.py` to verify and regenerate reporting. Large arrays/weights/count caches are regenerable and ignored; root code and small evidence are versioned.','']
 (ROOT/'V27_RESULTS.md').write_text('\n'.join(parts),encoding='utf-8')
 (ROOT/'V27_SUBMISSION.md').write_text('# V27 artifact decision\n\nNo V27 file is ready to upload. '+('The fixed finalist failed the prospective H1 control/treatment checks; no reserve was opened and no new package was generated.' if decision else 'Required tests or release checks are incomplete.')+'\n\nV24 remains the best recorded published submission (entry `81Y2yX7utdAau9rr3jCj`, overall 0.14299785521745675). Preserve its existing artifact; no repeat upload is needed. See [V27_RESULTS.md](V27_RESULTS.md) for measured comparisons and [V27_STRATEGY.md](V27_STRATEGY.md) for reproduction.\n',encoding='utf-8')
 readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8-sig');start='<!-- v27-control-anchored:begin -->';end='<!-- v27-control-anchored:end -->'
 block=start+'\n\n## V27 control-anchored neural response tests (2026-10-15)\n\nExecuted stage-by-stage sampling diagnostics and a matched panel-output neural fit, with nine paired population arms spanning fallback locking, query-control quantile calibration and relative neural response changes. '+('The locked finalist failed prospective H1 checks; all proposed pilot arms received held-family treatment tests. No V27 package is released, the donor reserve remains unopened, and V24 stays best recorded at 0.14299785521745675.' if decision else 'Execution/release validation is still incomplete; no upload-ready artifact is claimed.')+'\n\n[Methods and reproduction](V27_STRATEGY.md), [measured results and limitations](V27_RESULTS.md), [artifact decision](V27_SUBMISSION.md). Root code and small evidence are versioned; regenerable large arrays stay outside Git.\n\n'+end
 if start in text:
  a=text.index(start);bb=text.index(end,a)+len(end);text=text[:a]+block+text[bb:]
 else:
  i=text.index('\n');text=text[:i+1]+'\n'+block+'\n\n'+text[i+1:]
 readme.write_text(text,encoding='utf-8')
 try:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,axs=plt.subplots(1,2,figsize=(12,5));labels=sampling.stage.str.replace('_',' ');colors=['#246b45' if c==0 else '#a84936' for c in sampling.null_DE_calls];axs[0].barh(labels,sampling.null_DE_calls,color=colors);axs[0].invert_yaxis();axs[0].set_xlabel('Significant genes versus real NTC reference');axs[0].set_title('H1 control sampler ablation')
  if (A/'null_H1_scores.csv').exists():
   nf=pd.read_csv(A/'null_H1_scores.csv');means=nf.groupby('candidate').null_DE_calls.mean().reindex(ARMS);axs[1].barh(ARMS,means,color='#47789f');axs[1].invert_yaxis();axs[1].set_xlabel('Mean significant genes across three splits');axs[1].set_title('V27 matched control populations')
  fig.tight_layout();fig.savefig(A/'control_calibration.png',dpi=160);plt.close(fig)
 except ImportError:pass
 write(A/'report_manifest.json',{'plan_sha256':b.sha256(A/'plan.json'),'documents':{n:b.sha256(ROOT/n) for n in ['V27_STRATEGY.md','V27_RESULTS.md','V27_SUBMISSION.md']},'reporter_sha256':b.sha256(Path(__file__)),'completed_families':completed})
 print('V27 results documented',flush=True)
if __name__=='__main__':main()
