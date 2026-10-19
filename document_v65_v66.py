"""Document complete fixed V65/V66 experiments and the audited selected upload file."""
from pathlib import Path
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'experiments/v65_closed_response/analysis'
FOLDERS={65:'v65_closed_response',66:'v66_assay_strength'}

def main():
 decision=read(BASE/'release_selection.json');selected=decision['selected_version'];broad=read(BASE/'broader_finalist_diagnostic.json')
 for version,folder in FOLDERS.items():
  a=ROOT/'experiments'/folder/'analysis';mean=read(a/'mean_results.json');count=read(a/'count_assessment.json');audit=read(a/'results_audit.json');sensitivity=read(a/'paired_profile_sensitivity.json');assert audit['passed'];arm=count['selected']
  table=['| Response candidate | Profile error / V63 | LFC MSE / V63 | LFC MAE / V63 | PDS change | Eligible? |','| --- | ---: | ---: | ---: | ---: | :---: |']
  for name,s in mean['summary'].items():table.append(f"| {name} | {s['profile_error_ratio']:.9f} | {s['lfc_square_ratio']:.9f} | {s['lfc_absolute_ratio']:.9f} | {s['pds_change']:+.9f} | {s['eligible']} |")
  ct=['| Count arm | PDS | Fidelity | Reach | Jaccard | NMAE error | Profile squared error |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
  for name,s in count['summary'].items():ct.append(f"| {name} | {s['pds']:.9f} | {s['fidelity']:.9f} | {s['reach']:.9f} | {s['jaccard']:.9f} | {s['nmae']:.9f} | {s['profile_square_error']:.9f} |")
  ratio=count['paired_equal_study_profile_error_ratio'];bb=broad['summary']['v'+str(version)];failed=[k for k,v in count['checks'].items() if not v]
  st=['| Public study | Profile error / V63 | PDS change | Fidelity change | Reach change |','| --- | ---: | ---: | ---: | ---: |']
  for study in ('H1','Orion','CD4'):
   base=count['studies']['v63'][study];new=count['studies'][arm][study]
   st.append(f"| {study} | {new['profile_square_error']/base['profile_square_error']:.9f} | {new['pds']-base['pds']:+.9f} | {new['fidelity']-base['fidelity']:+.9f} | {new['reach']-base['reach']:+.9f} |")
  text=f'''# V{version} measured response and count results

**Status: {'selected for an explicitly experimental upload test' if version==selected and decision['experimental'] else 'selected by public guards, private score unknown' if version==selected else 'tested alternative; no production upload file generated'}.** V63 remains the best published project submission, overall 0.150422795. No V65/V66 leaderboard result is known.

## Fixed response screen

{chr(10).join(table)}

These are desired-mean diagnostics over the canonical first case for each of 30 unique target/context pairs. Contexts are averaged within study, then three studies equally. Profile ratios compare summed target errors within context to V63, excluding each intervention's own gene. PDS excludes all 300 panel targets. LFC errors have the inherited expression/support masks. Lower errors are better. These are repeatedly reused public datasets; no new independent biological confirmation is claimed.

Locked count candidate: `{arm}`. No alternate arm was promoted after count scores. V65 tests available-axis reclosure/perturbation closure; V66 is the strength-matched original-assay control. The normalized half-strength variants were rejected. A stronger correction in the original representation must be compared before crediting normalization for a gain. See [V65 protocol](V65_STRATEGY.md) and [V66 protocol](V66_STRATEGY.md).

## Complete paired integer-count tests

{chr(10).join(ct)}

{chr(10).join(st)}

These are external public study groups; they do not identify anonymous contexts A, B or C.

There are 240 comparison records: 30 unique target/context pairs, two technical generation seeds, four arms (V24, V63, primary and shuffled). Reused baseline files are byte-checked computations, not extra biology. Canonical target selection avoids repeated-label PDS competitors. Both generations use the same measured truth, so they are not independent replicates.

Public guard result: **{count['release_supported']}**. Failed checks: **{', '.join(failed) if failed else 'none'}**. Lower NMAE and profile errors are favorable; PDS/FID/REACH/JAC increase when favorable. HCT116's unresolved NMAE remains missing; the summary averages available studies for that metric. No official hidden anchors are estimated. The table's profile column averages absolute squared errors; the following paired ratio instead averages within-study candidate/baseline ratios, so it is not the ratio of those two table averages.

Paired equal-study generated-profile error / V63: **{ratio:.9f}** ({(ratio-1)*100:+.3f}%). Per-seed ratios: `{sensitivity['seed_ratios']}`. Leave-one-context-target-out range: **{sensitivity['loo_range'][0]:.9f} to {sensitivity['loo_range'][1]:.9f}**. Where the leave-one-out range crosses one, removing a single target/context reverses the aggregate gain. These are influence/technical-seed diagnostics, not confidence intervals or independent biological replication.

The serialized auditor independently replays fidelity, Jaccard, NMAE, call counts, PDS and profile error across all 240 records and 40 panels, and verifies reuse receipts and arrays. The legacy cache omits raw p values used for reach tie ordering; no full serialized reach replay is claimed. The selected release receives ten fresh prespecified first-target count/DE replays across all contexts/seeds, including an independent raw-p reach ordering check. Production integrity is checked separately.

## Measured count trade-off

![Study-specific profiles and raw-metric changes](experiments/v65_closed_response/analysis/count_tradeoff.png)

Each relative raw-metric change has its own scale. The percent bars must not be averaged or interpreted as official anchored-score contributions.

## Broader identity and expression check

The same existing raw extractions support 655 target/context pairs with at least 30 treated cells. On desired means before count generation, equal-study PDS change versus V63 is **{bb['pds_change']:+.9f}** and profile error ratio is **{bb['profile_ratio']:.9f}**. Profile scoring excludes own/unmeasured genes; PDS excludes the whole target panel. This broader result is diagnostic, not another selection set; no strength was changed to improve it.

The broader expression ratios are slightly worse than V63 for both variants, while the small fixed response/count cohorts show some gains. This disagreement and the failed count guards limit the evidence. No substantial rank improvement is demonstrated. Any released file is a bounded experimental test of a fixed candidate, not a replacement proven superior on the challenge. The unchanged generator still has the previously documented control-null DE calibration problem.

## Reproduction and provenance

Root scripts and exact input/implementation hashes are recorded in plan.json and production_plan.json; the whole-study source exclusions and poisoning checks pass. V66 also replays V63 bitwise at its original correction strength in every context. D3 remains unopened. Source studies, outputs and controls have been examined repeatedly in earlier versions; per-study source exclusions prevent recipient response leakage but do not create fresh model-selection data.

[V65 protocol](V65_STRATEGY.md) · [V66 protocol](V66_STRATEGY.md) · [All published submissions](V65_SUBMISSION_HISTORY.md) · [Selected submission and limitations](V65_V66_SUBMISSION.md). Compact receipts are saved under experiments/{folder}/analysis. The selected package receives official prep plus full payload/axis/provenance checks and 12 full-group production replays. These are integrity checks, not hidden accuracy tests. No automatic upload or Git push.
'''
  (ROOT/f'V{version}_RESULTS.md').write_text(text,encoding='utf-8')
 assert selected is not None,'No candidate qualified for an upload artifact'
 a=ROOT/'experiments'/FOLDERS[selected]/'analysis';ready=read(a/'submission_ready.json');audit=read(a/'independent_validation.json');assert audit['passed'] and audit['package_sha256']==ready['sha256'];assert read(a/'count_replay_audit.json')['passed'];assert read(a/'release_complete.json')['independent_integrity_passed'];candidate=decision['selected'];count=read(a/'count_assessment.json');primary=count['summary'][candidate['candidate']];baseline=count['summary']['v63']
 name='Closed-source assay response v65 experimental' if selected==65 else 'Assay-strength response v66 experimental'
 command=f'.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\{FOLDERS[selected]}\\prediction.vcc -m "{name}" --wait'
 report=f'''# Selected V{selected} submission: experimental test

Reviewed all 18 published submissions and tested four response-estimator variants across V65/V66. The fixed comparison selects **V{selected} `{candidate['candidate']}`** under the frozen rule. **V63 (0.150422795) remains our best published reference.** This new file has no private leaderboard result. Public all-guard support: **{candidate['all_public_guards_passed']}**. Failed checks: **{', '.join(candidate['failed_guards']) if candidate['failed_guards'] else 'none'}**.

## Why this file

The trajectory favors pooled intervention evidence, calibrated coverage, context-matched donor responses and the recent assay-representation correction; ordinary blending/stacking did not improve on V24. V65/V66 isolate source normalization from correction strength while retaining the same count generator. The selection rule was committed before V66 count outcomes: prefer all-guard passers, otherwise only bounded diagnostic candidates with at least two raw improvements and no study profile degradation over 3%, then the smallest paired generated-profile error ratio. See release_selection.json for every finalist's outcome.

Selected count-profile error / V63: **{candidate['profile_error_ratio']:.9f}**. Raw metric changes: PDS {primary['pds']-baseline['pds']:+.9f}; fidelity {primary['fidelity']-baseline['fidelity']:+.9f}; reach {primary['reach']-baseline['reach']:+.9f}; overlap {primary['jaccard']-baseline['jaccard']:+.9f}; NMAE error {primary['nmae']-baseline['nmae']:+.9f}. The generated-profile gain is concentrated in the held Orion study (HCT116/HEK293T); H1 and CD4 profile errors are slightly worse. These small reused-public-data trade-offs cannot forecast the official scaled overall. The 655-pair broader mean check is mixed and does not establish a clear improvement. Uploading measures this uncertainty; it is not a promised rank gain.

## Upload from Windows Command Prompt

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
{command}
```

Use this CMD syntax in the same terminal as previous submissions. If the wait is interrupted after upload, check the returned entry ID with `.\\.competition_env\\Scripts\\python.exe -m vcc status ENTRY_ID --wait`; do not upload a duplicate. This release has not been uploaded automatically.

## Artifact verification

- File: `experiments/{FOLDERS[selected]}/prediction.vcc`, **{ready['bytes']:,} bytes**.
- SHA-256: `{ready['sha256']}`.
- Official local prep passes: 360,000 cells, 18,533 genes, 900 context/target groups, 400 cells each.
- Every packaged sparse count/index/pointer and all axes match the generated source.
- All original control identities, 400 distinct templates/group, guides and exact libraries checked; all counts are nonnegative integers within cell limits.
- All 900 production response vectors replay with every prior-margin constraint checked; 12 full production groups regenerate exactly.
- Ten additional primary public populations freshly replay profile/DE arrays and independent reach tie ordering; complete serialized audits cover 480 comparison records across both finalists.

These checks establish file/provenance integrity, not predictive accuracy. Scientific source and compact receipts are locally committed; large data/arrays/packages are intentionally ignored by Git. No Git push.

[All 18 published submissions](V65_SUBMISSION_HISTORY.md) · [V65 results](V65_RESULTS.md) · [V66 results](V66_RESULTS.md) · [V65 protocol](V65_STRATEGY.md) · [V66 protocol](V66_STRATEGY.md).
'''
 (ROOT/'V65_V66_SUBMISSION.md').write_text(report,encoding='utf-8')
 readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8-sig');begin,end='<!-- v65-v66-release:begin -->','<!-- v65-v66-release:end -->'
 block=f'''{begin}
## V65/V66: submission history reviewed and new experimental test file

Reviewed **all 18 published submissions**; tested source normalization and a strength-matched assay correction with **480 paired count-comparison records** (shared baselines and repeated biology explicitly audited). Selected **V{selected} `{candidate['candidate']}`**, paired generated-profile error / V63 **{candidate['profile_error_ratio']:.9f}**. Public strict guards pass: **{candidate['all_public_guards_passed']}**; failures: {', '.join(candidate['failed_guards']) if candidate['failed_guards'] else 'none'}. The broader 655-pair desired-mean check is mixed. **V63 remains best published (0.150422795); this is an experimental upload test with unknown private accuracy.** Official packaging, full payload/provenance checks and deterministic replays pass. Nothing uploaded automatically.

[18-submission history and plot](V65_SUBMISSION_HISTORY.md) · [V65 results](V65_RESULTS.md) · [V66 results](V66_RESULTS.md) · [Selected file and CMD command](V65_V66_SUBMISSION.md). Older blocks below describe earlier decisions.
{end}
'''
 if begin in text:
  start=text.index(begin);finish=text.index(end,start)+len(end);text=text[:start]+block.rstrip()+text[finish:]
 else:text=block+'\n'+text
 readme.write_text(text,encoding='utf-8');print('Documented complete V65/V66 results and audited selected V'+str(selected)+' upload file.',flush=True)
if __name__=='__main__':main()
