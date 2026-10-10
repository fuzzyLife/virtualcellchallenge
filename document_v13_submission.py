"""Document only a completed, independently verified V13 package."""
import json
from pathlib import Path
import pandas as pd
from v13_methods import ROOT,A,OUT
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
def table(frame,columns):
 lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
 for _,r in frame.iterrows():
  lines.append('| '+' | '.join(f'{r[c]:.6f}' if isinstance(r[c],float) else str(r[c]) for c in columns)+' |')
 return '\n'.join(lines)
def main():
 ready=json.loads((A/'submission_ready.json').read_text(encoding='utf-8'));verified=json.loads((A/'independent_validation.json').read_text(encoding='utf-8'));model=json.loads((A/'reference_model.json').read_text(encoding='utf-8'));choice=json.loads((A/'count_selection.json').read_text(encoding='utf-8'));confirmation=json.loads((A/'confirmation_result.json').read_text(encoding='utf-8'))
 assert verified['passed'] and ready['sha256']==verified['package_sha256']
 folder=A/choice['validation_folder'];candidate=choice['candidate'];tables=[];failures=[]
 uncertainty=json.loads((A/'selected_uncertainty.json').read_text(encoding='utf-8'));interval_rows=[]
 for stage,record in uncertainty['results'].items():
  for metric,limits in record['conditional_DE_bootstrap_95pct'].items():
   interval_rows.append({'stage':stage,'metric':metric,'delta':record['point_family_balanced_delta'][metric],'low':limits['low'],'high':limits['high'],'families':record['eligible_family_count_per_metric'][metric]})
 intervals=table(pd.DataFrame(interval_rows),['stage','metric','delta','low','high','families'])
 for stage in ['development','confirmation']:
  scores=pd.read_csv(folder/f'{stage}_summary.csv');scores=scores[scores.reference_mode.eq('full_independent')&scores.candidate.isin(['v11',candidate])]
  tables.append('### '+stage.capitalize()+' (full independent controls)\n\n'+table(scores,['candidate']+FIELDS))
  assessment=json.loads((folder/f'{stage}_assessment.json').read_text(encoding='utf-8'))
  r=assessment['results'][candidate]
  for mode,record in r['modes'].items():
   failures.append(stage+' / '+mode+': '+(', '.join(k for k,v in record['checks'].items() if not v) or 'all guards pass'))
 status='Passed the fixed public development and confirmation guards. Hidden leaderboard improvement remains unverified.' if confirmation['confirmed_improvement'] else '**Exploratory candidate: local tests did not confirm improvement over V11. V11 remains the recommended best-supported submission.**'
 anchored=f"Public development anchored overall: {choice['source_anchored_overall']:.6f} versus V11 {choice['baseline_source_anchored_overall']:.6f} (difference {choice['source_anchored_overall']-choice['baseline_source_anchored_overall']:+.6f}). The raw-guard policy selects expression-error improvement while limiting other losses; this diagnostic does not demonstrate a better overall score."
 display={'v13_methods':'Shared response program transport v13','v13_robust_methods':'Robust source consensus transport v13','v13_quality_methods':'Assay-aware knockdown transport v13'}[choice['methods_module']]
 evaluator={'v13_methods':'evaluate_v13_submission.py','v13_robust_methods':'evaluate_v13_robust.py','v13_quality_methods':'evaluate_v13_quality.py'}[choice['methods_module']]
 command=f'.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v13_response_programs\\prediction.vcc -m "{display}" --wait'
 text=fr"""# V13 submission and validation

{status}

A new package is ready at `experiments/v13_response_programs/prediction.vcc`. Official prep and independent full-payload/replay checks passed. It has not been uploaded.

- Selected model: `{candidate}` from `{choice['methods_module']}`.
- Package size: {ready['bytes']:,} bytes ({ready['bytes']/2**30:.2f} GiB).
- SHA-256: `{ready['sha256']}`.
- Cells/genes: 360,000 by 18,533; 900 context-target groups with 400 cells each.
- Nonzero counts: {ready['prediction_validation']['nnz']:,}.
- Published reference: V11, overall 0.137003974; no V13 leaderboard result exists.

## Research and tests

See [V13_STRATEGY.md](V13_STRATEGY.md) for the primary literature, complete response alternatives, data coverage, validation design, correction record, and reproduction. We tested six low-rank variants, three robust-consensus variants, and three assay-aware knockdown variants at summary level. Two alternatives from each of the three classes advanced independently to generated-count development. One candidate was frozen before reserved confirmation, then compared once with V11. Development comprises 58 target-seed cases from three source families; confirmation comprises 30 cases from HCT116 and HEK293T only, with no H1 confirmation. Reused source cells, technical strata and paired seeds are not independent biological replication.

{(chr(10)*2).join(tables)}

Lower NMAE/expression error is better; higher PDS/fidelity/reach/Jaccard is better. These are **unanchored public-source metrics**, not competition scores. Both reference modes are retained in CSV. Some sampled400 panels lack an eligible NMAE contrast; they are reported as missing, never zero. Full independent development panels are complete. Missing confirmation panels, if any, are listed in the assessment records and prevent a full primary-guard pass. Any displayed NMAE average then covers eligible panels only; it is not evidence that the complete confirmation passed. A reporting-only adapter preserves these missing values without altering frozen predictions, metrics, or selection.

Guard results:

{chr(10).join('- '+f for f in failures)}

Conditional target-bootstrap intervals (2,000 Bayesian draws; candidate minus V11):

{intervals}

These intervals condition on the available public source families and do not measure uncertainty over unseen ABC lineages. Missing NMAE panels are excluded explicitly; zero change on an eligible subset is not complete validation. PDS and expression error are nonlinear panel metrics, so their paired family ranges are recorded separately in `selected_uncertainty.json`.

The failed low-rank stage remains recorded and was not relabeled as an improvement. {anchored}

Official public-source anchored comparisons are separate diagnostics; their small denominators and different contexts prevent interpreting them as leaderboard forecasts. Adaptive development selection increases optimism; frozen count confirmation is the final local check.

## Integrity and reproduction

The held-family mutation tests passed for all three response classes. All V11 baseline profiles, jackknifes and DE metrics replay exactly against independently generated V12 caches. The production response array was independently refitted and compared exactly. A separate post-freeze audit records small float64/float32 pooling differences outside the 65 KD-corrected targets (maximum 8.345e-7 log2 units across 235 other targets); those differences are part of the tested candidate, so this is not a bitwise single-change ablation. The payload audit compares every count, CSR index/pointer and axis; verifies all 360,000 sampled control identities, guides and depths; and exactly regenerates 12 complete groups ({12*400*18533:,} gene counts). Counts are nonnegative raw integers, with the required genes/order and cell/group sizes. File integrity does not establish biological prediction accuracy.

After completed development and selection:

```cmd
.\.competition_env\Scripts\python.exe audit_v13_baseline_replay.py
.\.competition_env\Scripts\python.exe select_v13_submission.py
.\.competition_env\Scripts\python.exe evaluate_v13_confirmation.py --workers 4
.\.competition_env\Scripts\python.exe confirm_v13_submission.py
.\.competition_env\Scripts\python.exe fit_v13_submission.py
.\.competition_env\Scripts\python.exe build_v13_submission.py --preflight --workers 4
.\.competition_env\Scripts\python.exe build_v13_submission.py --workers 4
.\.competition_env\Scripts\python.exe validate_v13_submission.py
.\.competition_env\Scripts\python.exe summarize_v13_evidence.py
.\.competition_env\Scripts\python.exe document_v13_submission.py
```

Run in a fresh checkout with the documented source data, or preserve the already frozen plan/model/package rather than overwriting it. Large source/count caches and packages are intentionally excluded from Git; root code and JSON/CSV evidence are tracked.

## Upload from Windows Command Prompt

```cmd
cd /d L:\promec\Animesh\Download\Virtual
{command}
```

An upload is a leaderboard experiment. If interrupted after scoring starts, use `vcc status <entry-id> --wait` instead of creating another entry. V11's package and its published result remain intact.
"""
 (ROOT/'V13_SUBMISSION.md').write_text(text,encoding='utf-8')
 readme=ROOT/'README.md';s=readme.read_text(encoding='utf-8');marker='<!-- V13_RESPONSE_TRANSFER -->'
 s=s.replace('no V12 package was generated and confirmation remains unopened.', 'no V12 package was generated. Confirmation remained unopened at V12 close; V13 later opened its reserved generated-count comparison, documented below.')
 if marker in s:s=s[:s.index(marker)].rstrip()+'\n'
 s+='\n'+marker+'\n\n## V13 response-transfer experiments and upload artifact\n\n'
 s+=f'Created `experiments/v13_response_programs/prediction.vcc` ({ready["bytes"]/2**30:.2f} GiB), using `{candidate}` after separate response-program, robust-consensus and knockdown-quality tests. Official packaging and independent complete-payload/provenance/replay audits passed. '
 s+=('The fixed public development and count-confirmation guards passed; a hidden leaderboard gain is unverified. ' if confirmation['confirmed_improvement'] else 'This is an **exploratory** artifact: the fixed local guards did not confirm an improvement over V11. V11 remains the best published and best-supported submission. ')
 s+='See [V13_STRATEGY.md](V13_STRATEGY.md) and [V13_SUBMISSION.md](V13_SUBMISSION.md) for literature, all outcomes, limitations, reproduction, hash, and Command Prompt upload command. No upload was performed.\n'
 readme.write_text(s,encoding='utf-8')
 build.write_json(A/'documentation.json',{'passed':True,'package_sha256':ready['sha256'],'report_sha256':build.sha256(ROOT/'V13_SUBMISSION.md'),'strategy_sha256':build.sha256(ROOT/'V13_STRATEGY.md'),'script_sha256':build.sha256(Path(__file__))})
 print('Completed V13 report:',model['status'],flush=True)
if __name__=='__main__':main()
