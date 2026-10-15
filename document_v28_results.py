"""Report actual V28 execution scope and artifact decision; no metric retuning."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v28_methods import ROOT,A,FINALIST,load_plan,write
from document_v27_results import table
import build_submission as b

def main():
 plan=load_plan();decision=json.loads((A/'release_decision.json').read_text()) if (A/'release_decision.json').exists() else None;parts=['# V28 results: global transport without depth strata','',('**No V28 submission is released.** The sole candidate failed prospective checks; V24 remains best recorded at0.14299785521745675.' if decision else '**No upload-ready V28 package is claimed.** Further validation or release checks remain required.'),'',f'Frozen plan SHA256: `{b.sha256(A/"plan.json")}`. V27 methods/outcomes remain unchanged. This follow-up was selected from spent control/sampler evidence, not V28 treatment scores.','']
 rows=[];completed=[]
 for key in ['H1','HCT116','HEK293T','CD4']:
  path=A/f'null_{key}_scores.csv'
  if not path.exists():continue
  f=pd.read_csv(path);completed.append('null_'+key)
  for context in f.family.unique():
   for arm in ['v24',FINALIST]:
    v=f[(f.family==context)&(f.candidate==arm)];rows.append({'context':context,'arm':arm,'mean_DE_calls':float(v.null_DE_calls.mean()),'max_DE_calls':int(v.null_DE_calls.max()),'detection_MAE':float(v.detection_MAE.mean()),'log_variance_MAE':float(v.log_variance_MAE.mean()),'max_mean_L1':float(v.mean_reference_L1.max())})
 parts+=['## Real-control calibration','',table(pd.DataFrame(rows),['context','arm','mean_DE_calls','max_DE_calls','detection_MAE','log_variance_MAE','max_mean_L1']) if rows else 'No completed control comparison yet.','','The null splits are cell-disjoint within cases and use independent actual/reference/template/calibration-fit populations. Three seeds reuse biology and are not independent studies. The first H1 case is spent V26/V27 development biology. Every output retains V24 gene totals and template libraries; final CPM drift is measured. These control diagnostics do not prove treatment or leaderboard accuracy.','']
 for key in ['H1','HCT116','HEK293T','CD4']:
  rr=A/f'null_{key}_result.json'
  if rr.exists():
   r=json.loads(rr.read_text());parts+=[f'{key} null gate passed: **{bool(r["passed"])}**. Failed checks: '+(', '.join(k for k,v in r['checks'].items() if not v) or 'none')+'.','']
 tables=[]
 for key in ['H1','HCT116','HEK293T','CD4']:
  pp=A/f'treatment_{key}_panels.csv'
  if pp.exists():tables.append(pd.read_csv(pp));completed.append('treatment_'+key)
 if tables:
  f=pd.concat(tables,ignore_index=True);summ=[]
  for context in f.family.unique():
   base=f[(f.family==context)&(f.candidate=='v24')];row=f[(f.family==context)&(f.candidate==FINALIST)]
   for metric in ['pds','fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']:
    a=float(base[metric].mean());bb=float(row[metric].mean());summ.append({'context':context,'metric':metric,'V24':a,'global_anchor':bb,'difference':bb-a})
  parts+=['## Paired held-family treatments','',table(pd.DataFrame(summ),['context','metric','V24','global_anchor','difference']),'','These are raw public-development metrics, not host anchor-scaled scores. PDS uses pooled expression profiles and is mathematically protected by the exact V24 gene totals; NMAE is nearly fixed by the CPM contract. Distribution changes mainly affect significance, direction/coverage ordering and jackknife expression-noise correction. Tiny mean-rounding differences should not be interpreted as substantive mean-predictor improvements.','']
  h=f[f.family=='H1'].groupby('candidate',sort=False)[['fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']].mean().reset_index();parts+=['## H1 neural/calibration diagnostics','',table(h,['candidate','fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']),'','The global raw anchor is the sole eligible candidate. The relative, panel, calibrated and shuffled arms are matched H1 diagnostics only; they cannot replace a failed finalist. All models reuse the unchanged H1-excluded masked/panel fits.','']
  for key in ['H1','HCT116','HEK293T','CD4']:
   rr=A/f'treatment_{key}_result.json'
   if rr.exists():
    r=json.loads(rr.read_text());parts+=[f'{key} intermediate treatment guards passed: **{bool(r["passed"])}**. Failed checks: '+(', '.join(k for k,v in r['checks'].items() if not v) or 'none')+'.','']
 if (A/'development_result.json').exists():
  r=json.loads((A/'development_result.json').read_text());parts+=['## Full development decision','',f'Passed: **{bool(r["passed"])}**. Failed checks: '+(', '.join(k for k,v in r['checks'].items() if not v) or 'none')+'.','','```json',json.dumps(r['changes'],indent=2),'```','']
 if (A/'independent_audit.json').exists():
  audit=json.loads((A/'independent_audit.json').read_text());parts+=['## Independent audit','',f'Passed: **{bool(audit["passed"])}**. Checked{audit["paired_cases_checked"]} paired cases and{len(audit["full_axis_replays"])} full-axis replays. Real-null checks are recomputed independently; hashes, actual-cell labels, unique/disjoint identities, nonnegative integer axes, exact cell libraries/gene totals and CPM drifts are verified. The first H1 global-anchor replay matches the preexisting diagnostic count digest. H1 neural training is unchanged from the prior128-original-row/662-fitting-mean audit.','']
 if decision:parts+=['## Release decision','',decision['reason'],'','No new package, upload or reserve opening occurred. Execution stopped at the first failed frozen required stage; the scope tables list only completed comparisons. No old experiment was overwritten and no diagnostic arm was promoted.','']
 else:parts+=['## Remaining release requirements','','A development pass alone does not produce an upload-ready file. The locked candidate still needs the fresh author-MD5-verified D3Rest reserve with all CD4 response information excluded, then official packaging and full-payload/provenance/replay audits.','']
 parts+=['## Reproduction and limitations','','See [V28_STRATEGY.md](V28_STRATEGY.md) and [V27_RESULTS.md](V27_RESULTS.md). Root code and small evidence are versioned; raw/model/count caches are regenerable and ignored. Run the pinned scientific Python, then `audit_v28_experiment.py` and `document_v28_results.py`. Frozen numerical plans/outcomes refuse overwrite; verified caches can resume.','','The development biology is reused. HCT116 lacks NMAE and CD4 D1 treatment populations have200 actual cells. Guide/batch overlap is explicitly archived. The global anchor transports control variation and fixed predicted means; it does not learn identified target-specific stochastic dynamics. Better null calibration does not guarantee improved treatment prediction, source-to-query transfer or leaderboard rank.','']
 (ROOT/'V28_RESULTS.md').write_text('\n'.join(parts),encoding='utf-8');(ROOT/'V28_SUBMISSION.md').write_text('# V28 artifact decision\n\nNo V28 file is ready to upload. '+(decision['reason'] if decision else 'Required validation or release checks remain incomplete.')+'\n\nV24 remains the best recorded published submission: entry `81Y2yX7utdAau9rr3jCj`, overall0.14299785521745675. Preserve its existing artifact; no repeat upload is needed. See [V28_RESULTS.md](V28_RESULTS.md) for measured scope and [V28_STRATEGY.md](V28_STRATEGY.md) for methods.\n',encoding='utf-8')
 readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8-sig');start='<!-- v28-global-transport:begin -->';end='<!-- v28-global-transport:end -->';block=start+'\n\n## V28 global-only count transport (2026-10-15)\n\nSeparately frozen follow-up to the localized V27 depth-stratum defect. Global integer transport and bounded CPM swaps retain V24 means without an extra measurement draw or depth-stratum repair. '+('The sole eligible candidate failed prospective checks; no package is released and the donor reserve remains unopened. V24 stays best recorded at0.14299785521745675.' if decision else 'Execution/release validation remains incomplete; no upload-ready artifact is claimed.')+'\n\n[Strategy and commands](V28_STRATEGY.md), [measured results/limits](V28_RESULTS.md), [artifact decision](V28_SUBMISSION.md). Earlier numeric models and packages are preserved.\n\n'+end
 if start in text:
  i=text.index(start);j=text.index(end,i)+len(end);text=text[:i]+block+text[j:]
 else:
  i=text.index('\n');text=text[:i+1]+'\n'+block+'\n\n'+text[i+1:]
 readme.write_text(text,encoding='utf-8')
 if rows:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  rr=pd.DataFrame(rows);fig,ax=plt.subplots(figsize=(10,5));contexts=list(rr.context.unique());xx=np.arange(len(contexts));a=rr[rr.arm=='v24'].set_index('context').loc[contexts];bb=rr[rr.arm==FINALIST].set_index('context').loc[contexts];ax.bar(xx-.18,a.mean_DE_calls,.36,label='V24');ax.bar(xx+.18,bb.mean_DE_calls,.36,label='Global anchor');ax.set_xticks(xx,contexts,rotation=15);ax.set_ylabel('Mean significant genes in control-vs-control test');ax.set_title('Executed real-control calibration; three splits per context');ax.legend();fig.tight_layout();fig.savefig(A/'control_calibration.png',dpi=160);plt.close(fig)
 write(A/'report_manifest.json',{'plan_sha256':b.sha256(A/'plan.json'),'documents':{n:b.sha256(ROOT/n) for n in ['V28_STRATEGY.md','V28_RESULTS.md','V28_SUBMISSION.md']},'completed':completed,'reporter_sha256':b.sha256(Path(__file__))});print('V28 results documented',flush=True)
if __name__=='__main__':main()
