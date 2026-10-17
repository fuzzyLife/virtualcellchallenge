"""Document the completed fixed V55 diagnostic; no scientific file mutation."""
from pathlib import Path
import json
import run_v55_mean50 as r
ROOT=r.ROOT;A=r.A

def main():
 from summarize_v55_studies import main as summarize
 summarize()
 assessment=r.read(A/'assessment.json');audit=r.read(A/'independent_evaluation_audit.json');production=r.read(A/'production_effect_audit.json');ready=r.read(A/'submission_ready.json');delivery=r.read(A/'final_delivery_audit.json')
 assert audit['passed'] and production['passed'] and delivery['passed'] and ready['official_prep_passed']
 assert r.read(A/'challenge_controls_audit.json')['passed']
 assert r.read(A/'payload_roundtrip_audit.json')['passed']
 vals=assessment['summaries'];seed=assessment['by_seed']
 lines=['# V55: fixed half-blend diagnostic submission','',
 'This is the previously unsubmitted V37 mean50 response, produced in a new namespace. The user authorized a diagnostic after reviewing local-to-official validation mismatch. Original V37 failed its REACH guard and remains rejected in its original records. V24 remains the best published submission; V55 has no known leaderboard score.','',
 '## Completed paired tests','',
 '44 existing public populations,30 distinct context-target pairs,2 emission seeds,2 arms:176 predicted-population comparisons.40 original development comparisons were reused by checksums and136 new predictions generated. This is outcome-informed robustness testing on reused biology, not independent confirmation. D3 remains unopened. Rest D1 actual cells number200; other treatment populations400.','',
 '| Raw public metric | V24 | Fixed mean50 | Change |','|---|---:|---:|---:|']
 for k in r.METRICS:
  a,b=vals['mean50'][k],vals['v24'][k]
  lines.append(f'| {k} | {b:.9f} | {a:.9f} | {a-b:+.9f} |' if a is not None and b is not None else f'| {k} | unavailable | unavailable | unavailable |')
 lines+=['','Lower expression/NMAE errors are favorable; other raw metrics increase when favorable. These values are not anchored leaderboard scores. Equal targets within context, contexts within study and studies overall; repeated target populations are not independent competitors. Expression uses the existing panel jackknife calculation and retains its tiny-denominator limitation.','',
 '## Seed and case sensitivity','', '| Seed offset | FID change | Expression ratio relative change | REACH change |','|---|---:|---:|---:|']
 for s,z in seed.items():
  a,b=z['mean50'],z['v24'];lines.append(f"| {s} | {a['fidelity']-b['fidelity']:+.9f} | {a['common_expression_ratio']/b['common_expression_ratio']-1:+.4%} | {a['reach']-b['reach']:+.9f} |")
 lines+=['',f"Fidelity improves under every context-target removal: **{assessment['improvement_support']['fidelity_leave_target_out_positive']}**.",
 f"Expression improves under every context-target removal: **{audit['expression_improves_every_leave_target_out']}**.",'',
 'The separate diagnostic technical bounds pass: **'+str(assessment['diagnostic_release_allowed'])+'**. These bounds permit an informative upload; they do not establish improvement or retrospectively relax the original V37 rule. Original strict REACH failure is retained. The full seed and30-case influence values are in assessment.json and independent_evaluation_audit.json.','',
 '## What is packaged','',
 'For each context/target/output gene, effect = float32(V24 +0.5*(submitted V36 retrieval_contrast25 -V24)). This averages log2 effects, not integer matrices. No new network is trained. Unsupported/equal parent coordinates retain V24 exactly. Target-own effects remain excluded from these vectors and the unchanged generator applies its0.1 repression factor once.','',
 'The same V24 residual10 emitter, guide-stratified sampling, original control libraries, production seed20261011 and dual-moment calibration are used. The original V24 and V36 files remain unchanged. Synthetic cell IDs inherited from the validated generator use a v36 prefix; context/target labels and official axes are explicit.','',
 f"Upload file: experiments/v55_mean50_diagnostic/prediction.vcc; **{ready['bytes']:,} bytes**.",f"SHA-256: {ready['sha256']}",'',
 'Independent audits verify all528 saved-array DE metric values,20 profile-expression panels,30 target-removal expression calculations,5 unchanged null comparisons, the full production response blend and exact fallback. Preflight tests24 full-axis groups. Full output validates360,000 cells by18,533 genes,900 context-target groups with400 cells each. Every row library/provenance is checked. Official prep passes; full archive decompression compares every CSR value/index/pointer and gene/context/target label. Original control-null inflation persists; this experiment does not fix it.','',
 '## Submit from Windows Command Prompt','', '~~~bat','cd /d L:\\promec\\Animesh\\Download\\Virtual',
 '.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v55_mean50_diagnostic\\prediction.vcc -m "Fixed V24 V36 half-blend diagnostic v55" --wait','~~~','',
 'No automatic upload or Git push was performed. This is a diagnostic upload with uncertain score, not a recommended replacement for the V24 reference based on current evidence.','',
 'Reproduction: [V55_MEAN50_STRATEGY.md](V55_MEAN50_STRATEGY.md); root scientific driver run_v55_mean50.py, independent auditors audit_v55_evaluation.py/audit_v55_mean50.py, official packaging adapter package_v55_mean50.py. Source/array/receipt hashes are recorded in analysis/plan.json and the final delivery manifest. Large files remain local and regenerable.','']
 study=r.read(A/'study_summary.json')
 extra=['## Study-specific trade-offs','','| Study | FID change | REACH change | Expression ratio change |','|---|---:|---:|---:|']
 for x in study['studies']:
  z=x['delta'];extra.append(f"| {x['study']} | {z['fidelity']:+.9f} | {z['reach']:+.9f} | {z['common_expression_ratio']:+.9f} |")
 extra+=['','CD4 drives most of the fidelity gain; Orion fidelity declines. Expression improvement is concentrated in Orion and its unstable HCT116 ratio denominator. These results do not establish gains across biological studies; the diagnostic status is retained. study_summary.json is a post-outcome descriptive summary and does not alter predictions or gates.','']
 at=lines.index('## What is packaged');lines[at:at]=extra
 (ROOT/'V55_MEAN50_RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
 marker='<!-- v55-mean50:begin -->'
 p=ROOT/'README.md';body=p.read_text(encoding='utf-8-sig')
 assert marker not in body
 first=body.find('\n')+1
 note='\n'+marker+'\n\n## V55 fixed half-blend diagnostic\n\nRetested the previously unsubmitted V37 mean50 on all44 existing populations and two seeds (176 comparisons,30 distinct context-target pairs). Original V37 rejection is preserved; V24 remains best published. Completed independent source/score/production/full-payload checks. This is an exploratory diagnostic, with no claimed leaderboard improvement. [Results, uncertainty and submit command](V55_MEAN50_RESULTS.md), [protocol](V55_MEAN50_STRATEGY.md).\n\n<!-- v55-mean50:end -->\n'
 p.write_text(body[:first]+note+body[first:],encoding='utf-8',newline='\r\n')
 names=['study_summary.json','challenge_controls_audit.json','plan.json','assessment.json','independent_evaluation_audit.json','production_effect_audit.json','reference_model.json','selection.json','release_decision.json','preflight.json','prediction_validation.json','submission_ready.json','payload_roundtrip_audit.json','final_delivery_audit.json']
 r.write(A/'final_manifest.json',{'diagnostic_submission':True,'uploaded':False,'leaderboard_score':None,'original_V37_strict_guards_passed':False,'D3_opened':False,'package_sha256':ready['sha256'],'package_bytes':ready['bytes'],'receipts_sha256':{n:r.sha(A/n) for n in names},'documentation_sha256':{n:r.sha(ROOT/n) for n in ['V55_MEAN50_RESULTS.md','V55_MEAN50_STRATEGY.md']},'document_generator_sha256':r.sha(Path(__file__))})
 print('V55 results, command, README and final manifest written.',flush=True)
if __name__=='__main__':main()
