"""Render measured V59 results without inventing a private leaderboard forecast."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v59_state_composition import *

def table(df):
 cols=list(df.columns);lines=['| '+' | '.join(map(str,cols))+' |','|'+'|'.join(['---']*len(cols))+'|']
 for row in df.itertuples(index=False,name=None):
  cells=[(f'{v:.6f}' if isinstance(v,(float,np.floating)) and np.isfinite(v) else ('missing' if pd.isna(v) else str(v))) for v in row];lines.append('| '+' | '.join(cells)+' |')
 return '\n'.join(lines)

def main():
 response=read(A/'response_assessment.json');count=read(A/'count_assessment.json');oracle=pd.read_csv(A/'oracle_summary.csv');rf=pd.read_csv(A/'response_summary.csv');cf=pd.read_csv(A/'count_cases.csv');summary=cf.groupby(['source','arm'])[['fidelity','reach','jaccard','nmae','profile_MSE']].mean().reset_index()
 primary=summary[summary.arm.isin(['v24','joint50'])];null=pd.DataFrame([read(p) for p in sorted((A/'counts').glob('*_null.json'))])[['source','residual10_DE','real_control_DE']]
 fig,axs=plt.subplots(1,2,figsize=(12,4.8));r=rf.pivot(index='source',columns='arm',values='MSE');r.plot.bar(ax=axs[0],rot=30);axs[0].set_ylabel('Mean squared LFC error');axs[0].set_title('Held-study mean prediction');axs[0].legend(fontsize=8)
 o=oracle.set_index('source');o[['zero','composition','within','joint','direct_mean']].plot.bar(ax=axs[1],rot=30);axs[1].set_ylabel('Independent-view reconstruction MSE');axs[1].set_title('Diagnostic treatment-split reconstruction');axs[1].legend(fontsize=8)
 fig.tight_layout();fig.savefig(A/'diagnostics.png',dpi=160);plt.close(fig)
 decision='The fixed candidate passes the recorded public support checks; subsequent production/package audit is required.' if count['release_supported'] else 'The fixed candidate fails the recorded public support checks. V24 remains the preferred published submission; no improved V59 package is justified.'
 text=f'''# V59 state-composition experiment: measured results

{decision}

## Scope and integrity

All {response['cases']} eligible target-context response cases were evaluated at five arms. All {count['comparisons']} count comparisons completed, including the separate state-cell diagnostic. These are reused public development studies, not a private competition score or fresh biological validation. D3 remains unopened.

The original V59 protocol and source files are preserved. Before any state fit or treatment score existed, a recorded correction closed reconstructed CPM profiles to one million. All candidates and oracle controls use it. See analysis/pre_score_correction.json and run_v59_closed_stage.py. Raw rows, architectures, shrinkage, sources, arms, cohorts and release gates did not change. The corrected canonical controller is run_v59_release.py.

Each recipient study is excluded from stable-gene selection, PCA, KMeans and donor mechanism fitting. HCT116 and HEK293T are one Orion study; both CD4 states are one study. NTCs alone fit per-lane normalization. The recipient mechanism archive contains treatment diagnostics, but the deployable transfer reads only its control-state counts/means, basal profile and measured mask. The independent boundary audit makes all recipient treatment fields unavailable and requires bitwise-identical predictions. Unsupported joint coordinates retain V24 exactly. All 300 panel genes are excluded from state inputs.

## Independent treatment-split reconstruction

The oracle fits each of two disjoint treatment/control views and predicts the other, using its control population. Lanes/views average within target first. It is diagnostic treatment evidence, never a deployable input. Lower MSE is better.

{table(oracle)}

## Whole-study response prediction

Primary joint50 is a fixed 50% bounded correction toward the learned mechanism. Composition50 adjusts only the estimated population component; pooled50 removes query composition changes. Shuffled_joint50 uses cyclic target-deranged donor mechanisms. No diagnostic was promoted after outcomes. Baselines are the frozen V58 arrays: held-CD4 comparisons are unchanged, while H1/Orion comparator baselines additionally exclude CD4.

{table(rf[['study','source','arm','MSE','MAE']])}

Response gate: **{response['passed']}**. Equal-study primary/V24 MSE ratio: **{response['equal_study_mean_ratio']:.6f}**.

## Integer-cell comparison

All arms use the same independent real reference controls per context/seed. Five mean arms use common 400-cell templates and the unchanged residual10 emitter. The state_cells50 diagnostic mixes uniformly sampled V24-edited controls with source-mechanism-edited, state-weighted controls. It uses exact integer library conservation and sampling without replacement; its abundance inclusion probabilities are approximate, not asserted exact density-ratio draws. Missing NMAE remains missing. Expression MSE excludes only the directly targeted coordinate; nonpanel MSE is retained separately in count_cases.csv. Neither is an official anchor-scaled score.

{table(summary)}

## Zero-intervention controls

{table(null)}

The state-cell zero-intervention contract returns unedited real NTC rows. Null calibration alone does not establish response prediction. Treatment state assignment uses post-perturbation RNA, so the decomposition cannot identify survival, proliferation or individual counterfactual trajectories.

## Release decision

{decision}

```json
{json.dumps(count['checks'],indent=2)}
```

![V59 diagnostics](experiments/v59_state_composition/analysis/diagnostics.png)

## Reproduce

[Original prospective strategy](V59_STATE_COMPOSITION_STRATEGY.md) and the pre-score correction receipt preserve the protocol. Use the corrected controller:

```bat
set "SCIENCE_PY=C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe"
"%SCIENCE_PY%" -X utf8 run_v59_release.py response
"%SCIENCE_PY%" -X utf8 audit_v59_state_composition.py boundary
"%SCIENCE_PY%" -X utf8 run_v59_release.py counts
"%SCIENCE_PY%" -X utf8 audit_v59_state_composition.py final
"%SCIENCE_PY%" -X utf8 document_v59_results.py
```

Plan SHA256: `{sha(A/'plan.json')}`. Correction SHA256: `{sha(A/'pre_score_correction.json')}`. Large numerical caches are listed with checksums in analysis/array_manifest.json and remain regenerable local artifacts. No leaderboard upload or Git push is performed.
'''
 (ROOT/'V59_STATE_COMPOSITION_RESULTS.md').write_text(text,encoding='utf-8')
 readme=ROOT/'README.md';old=readme.read_text(encoding='utf-8-sig');start='<!-- v59-state-composition:begin -->';end='<!-- v59-state-composition:end -->';block=f'''{start}
## V59: state composition and within-state responses

Completed {response['cases']} public response cases and {count['comparisons']} count comparisons, with whole-study exclusion, shuffled/pooled controls and independent treatment-split reconstruction. Public release support: **{count['release_supported']}**. D3 remains unopened. {decision}

[Measured results](V59_STATE_COMPOSITION_RESULTS.md) · [Prospective protocol](V59_STATE_COMPOSITION_STRATEGY.md). Corrected runner: `run_v59_release.py`; original scientific bytes and the pre-score CPM closure correction are preserved.
{end}
'''
 if start in old:old=old[:old.index(start)]+old[old.index(end)+len(end):].lstrip('\n')
 readme.write_text(block+'\n'+old,encoding='utf-8')
 write(A/'documentation.json',dict(report_sha256=sha(ROOT/'V59_STATE_COMPOSITION_RESULTS.md'),plot_sha256=sha(A/'diagnostics.png'),writer_sha256=sha(Path(__file__)),plan_sha256=sha(A/'plan.json'),release_supported=count['release_supported']))
 print(decision,flush=True)
if __name__=='__main__':main()
