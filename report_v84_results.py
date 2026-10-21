"""Report the frozen V84 experiments, including failed gates and oracle limitations."""
import platform
from importlib.metadata import version
import numpy as np
import v84_intervention_operator as m
from v57_pie import read,write,sha

def number(value,digits=6):
 return 'unavailable' if value is None else f'{value:.{digits}f}'

def main():
 means=read(m.A/'mean_results.json');selection=read(m.A/'selection.json');held=read(m.A/'held_essential_results.json');donors={s:read(m.A/f'donor_assessment_{s}.json') for s in [0,620000]}
 audits=['operator_audit.json','donor_exclusion_audit.json','donor_audit_0.json','donor_audit_620000.json']
 for name in audits:assert read(m.A/name)['passed'],name
 selected=selection['selected'];allowed=bool(selected) and all(d['passes'].get(selected,False) for d in donors.values())
 source_baselines=[]
 for arm in ['program25','diagonal25']:
  for destination in m.SOURCES[1:]:
   r=next(r for r in held['records'] if r['arm']==arm and r['destination']==destination)
   source_baselines.append(dict(arm=arm,destination=destination,ratio_to_identity=r['MSE']/r['identity_MSE'],ratio_to_constant=r['MSE']/r['constant_MSE']))
 decision=dict(release_allowed=allowed,selected_on_development=selected,primary='program25',development_passed=selection['strict_development_passed'],donor_primary_passes={str(s):d['passes']['program25'] for s,d in donors.items()},new_upload_artifact_generated=False,source_results=source_baselines,existing_best_published_overall=.1508821875498449,hosted_score_forecast=None,audits={n:sha(m.A/n) for n in audits},script_sha256=sha(__file__))
 lines=['# V84 results: measured cross-context response operators','',
  f'**New submission qualification: {allowed}.** Development selected `{selected}`. Both fixed operators failed the original expression-improvement requirement. The primary was nevertheless evaluated on the separate donor as prespecified. No new upload package was generated; V66/V82 remains the best published project result, overall **0.1508821875498449**. No new hosted score is claimed.','',
  '[Frozen strategy and reproduction](V84_STRATEGY.md) | [Research and project gaps](POST_V82_RESEARCH.md) | [Coverage correction](V83_COVERAGE_AUDIT.md)','',
  '## Source-task results need strong baselines','',
  'Train on 513 shared essential interventions; evaluate 125 held interventions on 6,087 nonpanel outputs. No competition target is in that training panel. Lower ratios are better. Copying K562 is the identity baseline; the constant baseline is the destination fitting-set average response.','',
  '| Operator | Destination | Error / copied K562 | Error / constant response |',
  '| --- | --- | ---: | ---: |']
 for r in source_baselines:lines.append(f"| {r['arm']} | {r['destination']} | {number(r['ratio_to_identity'])} | {number(r['ratio_to_constant'])} |")
 lines+=['','The primary improves on copying K562 but loses to the constant-response comparator in all three destinations. The diagonal operator only narrowly improves on the constant comparator in Jurkat. These results do not establish recovery of transferable perturbation identity. The operators were not refitted on held targets.','',
  '## Reused public development screen','',
  'All 655 target/context pairs were included, with legacy recipient-study response evidence excluded. Macro averages give H1, Orion and CD4 equal study weight. These are reused public cases, not independent confirmation.','',
  '| Candidate | Profile error / V66 | PDS change | Pass every original development criterion |',
  '| --- | ---: | ---: | --- |']
 for a,r in means['summary'].items():lines.append(f"| {a} | {number(r['profile_ratio'],9)} | {number(r['pds_change'])} | {all(means['checks'][a].values()) if a in means['checks'] else 'negative control'} |")
 lines+=['','Require profile ratio <=0.995, worst-study ratio <=1.03, PDS loss <=0.005, and better profile error than the separately fitted shuffled model. No arm qualifies. The primary passes the latter three requirements but leaves expression error effectively unchanged. The diagonal arm also fails discrimination and the shuffled comparison.','',
  '## Project-unseen D3 donor: final emitted counts','',
  'The 14 targets were selected from metadata before scoring. The whole CD4 response family and merged donor teachers were excluded. Query inference used 5,006 NTCs; all comparisons used a disjoint 4,994-cell reference and matched 400-cell actual/template populations. Each technical seed is assessed separately. The donor belongs to a previously used study; the seeds may overlap cells and are not independent biological replicates.','',
  'Raw public metrics are shown below. Smaller profile error and NMAE are better; larger PDS, fidelity, reach and Jaccard are better. These are not hosted anchor-scaled scores or a predicted overall score.','',
  '| Seed offset | Arm | Profile error | PDS | Fidelity | Reach | Jaccard | NMAE |',
  '| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
 for seed,d in donors.items():
  for a,r in d['summary'].items():lines.append(f"| {seed} | {a} | {number(r['profile_MSE'])} | {number(r['pds'])} | {number(r['fidelity'])} | {number(r['reach'])} | {number(r['jaccard'])} | {number(r['nmae'])} |")
 lines+=['','| Seed | Primary profile error change | Primary PDS change | Failed original donor checks |','| ---: | ---: | ---: | --- |']
 for seed,d in donors.items():
  base=d['summary']['v66'];v=d['summary']['program25'];fail=[k for k,passed in d['checks']['program25'].items() if not passed]
  lines.append(f"| {seed} | {100*(v['profile_MSE']/base['profile_MSE']-1):+.4f}% | {v['pds']-base['pds']:+.6f} | {', '.join(fail) or 'none'} |")
  lines+=[]
 lines+=['','The comparator is the V66 algorithm recomputed with all CD4 responses excluded, not the submitted ABC array. Complete per-target records, NMAE availability lists and individual guard decisions remain in the two donor assessment JSON files. No target, strength or threshold was selected from donor outcomes.','']
 oracle_path=m.A/'oracle_assessment.json'
 if oracle_path.exists():
  oracle=read(oracle_path);audit=read(m.A/'oracle_audit.json');assert audit['passed'];decision['oracle_audit_sha256']=sha(m.A/'oracle_audit.json')
  lines+=['## Conditional diagnostic: mean versus population shape','',
   'After the fixed primary failed, the presaved protocol selected ADNP, MAP2K2, MYRF, PCGF1, USP28 and ZNF25. Each target uses disjoint 200-cell treatment fitting/scoring halves. Treatment labels enter only these non-deployable diagnostic populations.','',
   '| Diagnostic population | Profile error | Fidelity | Reach | Jaccard | NMAE | CPM contract holds for all six |',
   '| --- | ---: | ---: | ---: | ---: | ---: | --- |']
  for a,r in oracle['summary'].items():lines.append(f"| {a} | {number(r['profile_MSE'])} | {number(r['fidelity'])} | {number(r['reach'])} | {number(r['jaccard'])} | {number(r['nmae'])} | {r['CPM_contract_valid']} |")
  lines+=['','Within each mean pair, gene totals, row totals and raw profiles match exactly: `baseline` versus `predicted_mean_fitted_shape`, and `real_fit` versus `fitted_mean_old_shape`. The independent audit replays profiles, BH and direction metrics across all 24 populations and checks 192 original tied-rank coordinates. A failed CPM contract means that a projected population cannot be treated as isolating shape on every statistic.','',
   'Across mean pairs, fitted treatment depth and projection behavior also differ. Measured fitting-half means are noisy; six targets and 200-cell scoring halves do not establish a universal bottleneck. Oracle populations cannot train or select an upload model.','']
 if oracle_path.exists():
  ss=oracle['summary'];base=ss['baseline'];mean=ss['fitted_mean_old_shape'];real=ss['real_fit']
  lines+=['## What this changes in our priorities','',
   f"On this limited diagnostic, supplying measured fitting-half means reduces profile error by {100*(1-mean['profile_MSE']/base['profile_MSE']):.2f}% and raises raw fidelity from {base['fidelity']:.6f} to {mean['fidelity']:.6f}. That is evidence of response-mean headroom in these six cases, not evidence that our available inference inputs can recover those means.",'',
   f"DE overlap still rises from {mean['jaccard']:.6f} with fitted means and the old shape to {real['jaccard']:.6f} using real fitting-half cells. Matching means therefore leaves distributional errors. The fitted-shape-at-baseline-means branch fails its CPM contract, so its poor scores cannot isolate a pure shape effect.",'',
   'The next justified modeling target is accurate target-specific expression change in the scoring statistics, demonstrated on a whole held-out study and against constant/shuffled controls before adding a new population generator. More source data must first pass explicit target coverage or paired-response transfer tests. Another blend or dispersion setting without that evidence is unlikely to resolve the observed gap.','']
 lines+=['## Integrity and preserved failures','',
  '- Both donor seeds have 42 comparison populations, 84 total, with 672 independently checked tied-rank coordinates. Independent replay verifies BH, profiles, panel PDS and DE direction metrics.',
  '- The exclusion audit checks 28 lane-matched views and bitwise invariance after poisoning all excluded CD4 responses, errors, doses, masks and source profiles. Access to the CD4 teacher is forbidden in the replay.',
  '- Independent operator checks verify the normal equations, fitting split, scales, diagonal formula, held-target errors, zero-input identity and frozen source/code hashes.',
  '- V83 direct augmentation stopped at zero target coverage. The earlier V17 finding and initial parser engineering corrections are retained. Missing targets were never manufactured as biological zero votes.',
  '- The scDEBART code audit checks reviewed attention classes and configuration at a pinned commit. No released checkpoint was executed, and no claim is made about the code used for its published benchmarks.',
  '- Large source files, arrays, count populations and logs are ignored by Git. Compact provenance, code, fixed plans, per-case scores and audits are tracked. Original scientific releases were not rewritten.','',
  '## Reproduce and interpret','',
  'Use the commands in `V84_STRATEGY.md`, then run `report_v84_results.py`. The execution protocol records the model/code/input hashes before recipient evaluation. `finish_v84_donor.py` resumes completed case files only after checking their array hashes.','',
  'The new paired-evidence approach was implemented and tested rather than assumed to work from a paper. It does not currently justify replacing V66. See `POST_V82_RESEARCH.md` for the evidence gap and the distinction between our task and published benchmarks.']
 (m.ROOT/'V84_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 write(m.A/'release_decision.json',decision)
 write(m.A/'runtime_versions.json',dict(python=platform.python_version(),packages={p:version(p) for p in ['numpy','scipy','pandas','scikit-learn','h5py','pyarrow','anndata','torch','threadpoolctl','requests']},platform=platform.platform()))
 print('V84 final qualification',allowed,'report written',flush=True)
if __name__=='__main__':main()
