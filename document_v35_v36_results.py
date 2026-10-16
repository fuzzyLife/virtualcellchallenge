"""Document completed V35/V36 scientific tests and the requested upload artifact."""
import json
from pathlib import Path
import v35_methods as m

def read(p):return json.loads(p.read_text())
def main():
 a=read(m.A/'summary_results.json');b=read(m.P36/'summary_results.json');means=read(m.P36/'mean_results.json');selection=read(m.P36/'selection.json');ready=read(m.P36/'submission_ready.json')
 for folder in [m.A,m.P36]:assert read(folder/'independent_audit.json')['passed']
 assert read(m.P36/'independent_count_audit.json')['passed'] and ready['official_prep_passed']
 assert read(m.P36/'hybrid/independent_audit.json')['passed']
 assert read(m.P36/'contrast/independent_audit.json')['passed']
 assert read(m.P36/'production_effect_audit.json')['passed']
 names={'v35':'Experimental retrieval residual v35','v36':'Population-view response mean v36','v36_dispersion':'Population-view mean and dispersion v36','v35_population_dispersion':'Retrieval mean with population dispersion v36','target_centered25':'Target-centered experimental retrieval v36','retrieval_contrast25':'Retrieval-contrast response correction v36'}
 submission_name=names[selection['selected_arm']]
 bs=chr(92);directory=bs.join(['L:','promec','Animesh','Download','Virtual']);exe=bs.join(['.','.competition_env','Scripts','python.exe']);prediction=bs.join(['.','experiments','v36_population_views','prediction.vcc']);upload_command='cd /d '+directory+'\n'+exe+' -m vcc submit '+prediction+' -m \"'+submission_name+'\" --wait'
 selected_mean=selection['raw_mean_summaries'].get(selection['selected_arm'],selection['raw_mean_summaries']['v35' if selection['model_version']==35 else 'v36'])
 assert read(m.P36/'final_delivery_audit.json')['passed']
 rows=[]
 for v,r in [(35,a),(36,b)]:
  s=r['summary']['learned'];rows.append(f"| V{v} | {(s['mse_ratio']-1)*100:+.3f}% | {(s['mae_ratio']-1)*100:+.3f}% | {'Pass' if r['passes'] else 'Fail'} |")
 raw=[]
 for v in ['v35','v36']+[a for a in ['target_centered25','retrieval_contrast25'] if a in selection['raw_mean_summaries']]:
  r=selection['raw_mean_summaries'][v];raw.append(f"| {v.upper()} | {(r['mse_ratio']-1)*100:+.3f}% | {(r['mae_ratio']-1)*100:+.3f}% | {(r['cpm_mse_ratio']-1)*100:+.3f}% |")
 counts=[];base=selection['count_summaries']['v24']
 for arm in [a for a in ['v35','v36','v36_dispersion','v35_population_dispersion','target_centered25','retrieval_contrast25'] if a in selection['count_summaries']]:
  r=selection['count_summaries'][arm];fmt=lambda k: 'unresolved' if r[k] is None or base[k] is None else f"{r[k]-base[k]:+.6f}"
  err='unresolved' if r['common_expression_ratio'] is None else f"{(r['common_expression_ratio']/base['common_expression_ratio']-1)*100:+.3f}%"
  counts.append(f"| {arm} | {fmt('pds')} | {fmt('fidelity')} | {fmt('reach')} | {fmt('jaccard')} | {fmt('nmae')} | {err} |")
 null=read(m.P36/'null_results.json');nt=[]
 for c in sorted({r['family'] for r in null['records']}):
  rr={r['arm']:r for r in null['records'] if r['family']==c};nt.append(f"| {c} | {rr['v24']['real_null_DE_calls']} | {rr['v24']['null_DE_calls']} | {rr['v35']['null_DE_calls']} | {rr['v36']['null_DE_calls']} | {rr['v36_dispersion']['null_DE_calls']} |")
 coverage=[]
 for outer in m.OUTERS:
  r=read(m.A/'models'/(outer+'.json'));v=read(m.P36/'models'/(outer+'.json'));coverage.append(f"| {outer} | {r['outputs']} | {v['extra']['population_groups']} |")
 text='''# V35–V36 results and submission

V24 remains the published reference at overall 0.14299785521745675. These are reused public development tests, not private anchored leaderboard scores. A new experimental package was explicitly requested and prepared; it is not a confirmed better model.

## Gene-level summary response tests

| Model | Squared LFC error change vs V24 | Absolute LFC error change | Strict summary gate |
|---|---:|---:|---|
'''+ '\n'.join(rows)+'''

Whole H1, whole Orion (HCT116+HEK293T), and whole CD4 exclusions are averaged equally. The learned retrieval is compared with response identities shuffled within each source and fixed attention at inference. These controls are ablations, not separately trained negative-control models. V35 learns detailed experimental retrieval and a bounded residual; V36 adds independent raw-population views and a dispersion head. Neither is an exact reproduction of the cited papers.

## Coverage

| Excluded family | Shared modeled output genes | Eligible raw population groups, fitting and validation |
|---|---:|---:|
'''+ '\n'.join(coverage)+'''

Population groups repeat interventions across contexts; repeated sampled views are not independent biological experiments. In particular, the small Orion-excluded population branch limits evidence for broad transfer. Other genes retain V24 transfer-effect values exactly; normalization and sampling can change their emitted counts.

## Actual-cell mean tests

| Model | Squared LFC error change | Absolute LFC error change | CPM mean squared error change |
|---|---:|---:|---:|
'''+ '\n'.join(raw)+'''

All 44 previously spent raw mean cases were evaluated, with the matching family's responses excluded. Different denominators mean this table and the summary table are not interchangeable. No fresh held-out lineage or D3 donor result is claimed.

## Raw count diagnostics

| Arm | PDS change | FID change | REACH change | JAC change | NMAE change | Relative expression error change |
|---|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(counts)+'''

PDS/FID/REACH/JAC increases are favorable; NMAE/expression-error decreases are favorable. Twenty treatment cases: first four cases from each of five fixed public contexts. D1 treatment cases use 200 cells and are lower-power diagnostics; other cases use 400. V35/V36 means first use unchanged V24 count emission; the final arm separately tests learned gene-wise dispersion weights. Undefined values are never replaced with zero. Metrics average defined contexts within each study; unresolved studies fail their guard.

| Context | Actual null DE calls | V24 | V35 | V36 | V36 dispersion |
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(nt)+f'''

Null fidelity is not treatment recovery. Existing null artifacts are explicitly retained in the report; a plausible-looking cell cloud is not sufficient validation.

## Selected artifact

- Selected arm: **{selection['selected_arm']}**.
- Status: **{selection['status']}**.
- Selected mean error changes versus V24: squared LFC {(selected_mean['mse_ratio']-1)*100:+.3f}%, absolute LFC {(selected_mean['mae_ratio']-1)*100:+.3f}%, CPM {(selected_mean['cpm_mse_ratio']-1)*100:+.3f}%.
- Strict public improvement guards: **{selection['strict_public_guards_passed']}**.
- New blind confirmation: **False**. D3 reserve opened: **False**.
- Package: `experiments/v36_population_views/prediction.vcc`.
- Bytes: {ready['bytes']:,}; SHA-256: `{ready['sha256']}`.
- Official prep/validation: passed. Automatically uploaded: False.
- The original V24 package and champion record remain unchanged.

The secondary V35-mean/V36-dispersion component check was added after summary, raw-mean and partial count outcomes were visible. It uses identical cases, fixed weights and a separately frozen rule, and remains explicitly exploratory. Original primary selection is preserved in `primary_selection.json`.

A second exploratory correction was added after the full primary outcomes exposed a large PDS/expression regression. Two fixed quarter-strength transforms were tested: remove the target-shared correction, or retain genuine-minus-shuffled retrieval dependence. Complete 300 predictions define centering; no treatment labels enter the transforms. The same 44 means and 20 count cases were rerun for each. `precontrast_selection.json` preserves the prior decision; contrast results and audits remain separate.

The pre-outcome primary rule chooses between the two experimental mean candidates using equal-study raw mean error, preferring count noninferiority when available. Dispersion is adopted only if its separate predefined comparison passes. If strict gains are absent, the package remains exploratory; no leaderboard improvement is asserted.

## Verification and reproduction

Seven boundary tests passed. All 6 outer models' source/family/target exclusions and saved summary scores were audited. A separately written NumPy implementation checks the neural forward. Explicit disjoint cell identities are checked for each source/target's first pair; all pairs are asserted during construction. All 100 primary count-arm DE-set summaries and 20 primary treatment PDS/expression panels were replayed independently. The exploratory combination and contrast arms add 60 DE-set summaries and 15 profile panels, also independently replayed from saved arrays. Reach's p-value tie order and complete raw Wilcoxon computation were inherited from the frozen scorer, not independently regenerated.

Production preflight caught and repaired object-string label serialization before count generation. Initial files and generator source are preserved; numeric predictions were verified exactly unchanged without unpickling. The corrected adapter and independent production effect audit are recorded. Full count/library/gene/order/density validation and official lossless packaging are recorded in analysis receipts. A separate final audit independently checks every one of the 360,000 row libraries and its sampling-provenance/obs alignment. Large model/data/prediction arrays are ignored in Git; root scripts and metadata document their regeneration. See [strategy and exact commands](V35_V36_STRATEGY.md), [safe numerical replay](reproduce_v35_v36.py) and `requirements-v35-science.txt`.

Original plans and code are preserved. A separate pre-outcome selector amendment corrected the scorer field name from a nonexistent expression-error label to `common_expression_ratio`, without changing thresholds or the comparison. No numeric result was available when that initial field-name correction was made. A later separate contrast-selector repair handles undefined HCT116 NMAE using the unchanged primary aggregation; original predictions, failed selector and logs are preserved. One fixed training seed and one fixed emission seed per case were used; small gains are not statistical proof of generalization.

Windows Command Prompt upload command:

```bat
{upload_command}
```
'''
 published_path=m.P36/'leaderboard_result.json'
 published=None
 if published_path.exists():
  published=read(published_path);prior=read(m.ROOT/'experiments/v24_cd4_response/analysis/leaderboard_result.json')
  before=dict(overall=prior['overall'],**prior['metrics']);after=dict(overall=published['overall'],**published['metrics'])
  comparison='\n'.join(f"| {k.upper()} | {before[k]:.9f} | {after[k]:.9f} | {after[k]-before[k]:+.9f} |" for k in before)
  publication=f"""## Published leaderboard result

Entry `{published['entry_id']}` is **published**, with rank **{published['rank_at_user_report']}** at the user's report. Source: user-provided CLI output; no independent API confirmation is claimed. Same validation panel and anchors as V24.

| Scaled score (higher is better) | V24 | V36 | Change |
|---|---:|---:|---:|
{comparison}

V24 remains preferred. The overall decline is **{published['changes_from_v24']['overall']:+.9f}**. Expression accuracy and DE overlap improve, but PDS, NMAE, FID and especially REACH decline. Rank changes across dates also reflect other teams' submissions. This result tests the selected conservative retrieval contrast; the full population-view model was rejected locally and was excluded from this package. Small score differences from one submission do not establish statistical significance.

The generation receipts below retain their original pre-upload status and selection evidence. The published result is stored separately in `experiments/v36_population_views/analysis/leaderboard_result.json`; no training, model, selection or audit receipt was rewritten.

"""
  text=text.replace('## Gene-level summary response tests',publication+'## Gene-level summary response tests',1)
  text=text.replace('Windows Command Prompt upload command:', 'Original Windows Command Prompt upload command (entry already published; no repeat upload needed):')
 (m.ROOT/'V35_V36_RESULTS.md').write_text(text,encoding='utf-8')
 readme=m.ROOT/'README.md';s=readme.read_text(encoding='utf-8');start='<!-- v35-v36-model-tests:begin -->';end='<!-- v35-v36-model-tests:end -->';entry=f'''{start}

## V35/V36 retrieval and population-view tests (2026-10-16)

Implemented experimental gene-response retrieval with a bounded V24 residual, then independent raw-population-view mean/dispersion training. Both completed whole-family validation, 44 actual-cell mean cases, 20 count cases and 5 control-null cases. Summary squared LFC error changes: V35 {(a['summary']['learned']['mse_ratio']-1)*100:+.3f}%, V36 {(b['summary']['learned']['mse_ratio']-1)*100:+.3f}%. Negative error changes are favorable; these are reused public development results, not private VCC scores.

The requested upload package is [prediction.vcc](experiments/v36_population_views/prediction.vcc), selected arm `{selection['selected_arm']}`. **{selection['status']}**. Strict public guards passed: {selection['strict_public_guards_passed']}; new blind confirmation: False. V24 remains the published champion at 0.14299785521745675. Its package is unchanged; D3 remained sealed; no upload was performed.

Seven targeted tests and independent inference, score, split and count audits passed. Official packaging and a separate 360,000-row provenance/library audit also passed. [Results, limits and upload command](V35_V36_RESULTS.md), [strategy and commands](V35_V36_STRATEGY.md), [safe reproduction runner](reproduce_v35_v36.py). Scripts stay in the root; large arrays/packages remain outside Git.

{end}
'''
 if published is not None:
  entry=entry.replace('Its package is unchanged; D3 remained sealed; no upload was performed.', 'Its package is unchanged; D3 remained sealed. The user subsequently uploaded V36; entry `'+published['entry_id']+'` is published at rank '+str(published['rank_at_user_report'])+' (reported at submission), overall **'+f"{published['overall']:.9f}"+'**, change versus V24 **'+f"{published['changes_from_v24']['overall']:+.9f}"+'**. Expression accuracy and DE overlap improve, but the other four scaled metrics decline, especially REACH. V24 remains the preferred submission. Original pre-upload receipts remain unchanged; the published result is recorded separately.')
 if start in s:
  i=s.index(start);j=s.index(end,i)+len(end);s=s[:i]+entry+s[j:]
 else:s=s.replace('# Virtual Cell Challenge 2026\n','# Virtual Cell Challenge 2026\n\n'+entry,1)
 readme.write_text(s,encoding='utf-8');print('V35/V36 reports updated',flush=True)
if __name__=='__main__':main()