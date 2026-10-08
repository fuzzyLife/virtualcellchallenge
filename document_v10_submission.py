"""Render V10 facts only after official packaging and independent audit pass."""
import json
import pandas as pd
import build_submission as build
from v10_methods import ROOT,A,OUT

def table(frame,columns):
 result='| '+' | '.join(columns)+' |\n| '+' | '.join(['---']*len(columns))+' |\n'
 for _,row in frame.iterrows():
  values=[f'{row[c]:.6f}' if isinstance(row[c],float) else str(row[c]) for c in columns]
  result+='| '+' | '.join(values)+' |\n'
 return result

def main():
 ready=json.loads((A/'submission_ready.json').read_text(encoding='utf-8'))
 audit=json.loads((A/'independent_validation.json').read_text(encoding='utf-8'))
 model=json.loads((A/'reference_model.json').read_text(encoding='utf-8'))
 selection=json.loads((A/'count_selection.json').read_text(encoding='utf-8'))
 confirmation=json.loads((A/'confirmation_result.json').read_text(encoding='utf-8'))
 corrected=json.loads((A/'corrected_depth_validation.json').read_text(encoding='utf-8'))
 links=json.loads((A/'linked_provenance_validation.json').read_text(encoding='utf-8'))
 assert links['passed'] and links['package_sha256']==ready['sha256']
 depth=json.loads((A/'production_depth_audit.json').read_text(encoding='utf-8'))
 assert ready['official_prep_passed'] and audit['passed'] and depth['passed']
 assert depth['package_sha256']==ready['sha256'] and depth['verified_cells']==360000
 assert ready['sha256']==build.sha256(OUT/'prediction.vcc')==audit['package_sha256']
 assert selection['candidate']==model['candidate']==confirmation['candidate']
 output=f"# V10 verified submission\n\nSelected **{model['candidate']}**; response `{model['response']}`, emission `{model['generator']}`. This is an unsubmitted experimental prediction. V4 remains the best published submission until another official result beats it.\n\n"
 output+='See [V10_STRATEGY.md](V10_STRATEGY.md) for the diagnosis, public approaches, new response pooling, population transport and explicit neutral-control bias.\n\n'
 output+='## Corrected validation of the frozen candidate\n\n'
 output+='The original compact validation adapter indexed depths from the beginning of the control population rather than the recorded 400-cell template. Production indexes the full control matrix correctly. The separate root audit corrects the adapter, verifies every generated validation depth against its actual template, tests both seeds with 400 and full independent reference controls, and preserves the original snapshots. No model parameter or selection was changed. These are revisited public development targets, not new blind confirmation.\n\n'
 output+=table(pd.read_csv(A/'corrected_depth_summary.csv'),['reference_mode','candidate','pds','fidelity','reach','jaccard','nmae','common_expression_ratio'])+'\n'
 for mode,details in corrected['diagnostics'].items():
  difference=details['paired_six_metric_utility'][selection['candidate']]-details['paired_six_metric_utility']['v4_legacy']
  output+=f"Reference `{mode}`: original-style fidelity/PDS/expression guards passed **{details['guard_passed']}**; paired six-metric utility difference **{difference:.6f}**.\n\n"
 for stage in ['calibration','confirmation']:
  frame=pd.read_csv(A/f'count_{stage}_summary.csv')
  output+=f'## Original public {stage} snapshot (before depth-adapter repair)\n\n'+table(frame,['candidate','pds','fidelity','reach','jaccard','nmae','common_expression_ratio','utility'])+'\n'
 output+=f"The original historical confirmation snapshot passed guards: **{confirmation['guard_passed']}**; positive paired six-metric proxy improvement: **{confirmation['confirmed_improvement']}**; paired utility difference **{confirmation['utility_difference']:.6f}**. No candidate was selected or retuned using confirmation. H1 is reused development; HCT/HEK count targets are new, with public summaries previously available. These figures are not predictions of the hidden anchored score.\n\n"
 output+=f"Original snapshot full target-mean fidelity: **{confirmation['full_target_mean_fidelity']['v4_legacy']:.6f} -> {confirmation['full_target_mean_fidelity'][selection['candidate']]:.6f}**. On identical finite cases: **{confirmation['paired_finite_fidelity']['v4_legacy']:.6f} -> {confirmation['paired_finite_fidelity'][selection['candidate']]:.6f}**. The paired proxy uses {confirmation['common_six_metric_panels']} common complete panels; unresolved panels are not silently filled.\n\n"
 output+='![Corrected public validation comparisons](experiments/v10_population_transport/analysis/validation_comparison.png)\n\n'
 if (A/'corrected_depth_family_fidelity.csv').exists():
  output+='## Fidelity by source family (corrected adapter)\n\n'+table(pd.read_csv(A/'corrected_depth_family_fidelity.csv'),['reference_mode','family','unique_targets','baseline_fidelity','candidate_fidelity','difference'])+'\n'
  output+='The full-reference gain is uneven: H1 and HEK293T improve, while HCT116 declines. Reference size materially changes the apparent benefit. These public targets cannot establish transfer to the three hidden contexts.\n\n'
 output+='## Neutral-control diagnostic\n\n'+table(pd.read_csv(A/'decoder_neutral_summary.csv'),['generator','max_false_de_calls','max_false_call_fraction','mean_variance_ratio','neutral_guard'])+'\n'
 output+='Compact emission has substantial control-to-control DE bias and fails the former neutral-control preservation gate. The package is an experimental population model. Do not interpret synthetic dispersion as predicted biological variance.\n\n'
 output+=f"## Verified upload artifact\n\n- File: `experiments/v10_population_transport/prediction.vcc`\n- Bytes: {ready['bytes']:,}\n- Shape: 360,000 cells x 18,533 genes; 900 target/context groups x 400 cells.\n- Stored values: {ready['prediction_validation']['nnz']:,}.\n- Official prep, independent full-payload comparison and all 360,000 production-depth checks: passed.\n- Exactly replayed groups: {audit['exact_replayed_groups']}; replayed gene counts: {audit['exact_replayed_gene_counts']:,}.\n- SHA-256: `{ready['sha256']}`\n- Uploaded: false. Hidden leaderboard performance: unknown.\n\n"
 slash=chr(92)
 cwd=slash.join(['L:','promec','Animesh','Download','Virtual'])
 exe=slash.join(['.','.competition_env','Scripts','python.exe'])
 file=slash.join(['.','experiments','v10_population_transport','prediction.vcc'])
 output+=f'## Submit from Windows Command Prompt\n\n```bat\ncd /d {cwd}\n{exe} -m vcc submit {file} -m "Pooled intervention population transport v10" --wait\n```\n\n'
 output+='If Ctrl+C interrupts waiting after upload, use `vcc status <entry-id> --wait` to follow the existing scoring job instead of creating a duplicate.\n'
 (ROOT/'V10_SUBMISSION.md').write_text(output,encoding='utf-8')
 readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8');marker='<!-- V10_VERIFIED_SUBMISSION -->'
 std=corrected['diagnostics']['sampled400'];full=corrected['diagnostics']['full_independent']
 section=f"{marker}\n\n## V10 upload candidate\n\nV10 is ready at `experiments/v10_population_transport/prediction.vcc`. Selected `{model['candidate']}`; corrected sampled-control guard passed: {std['guard_passed']}; corrected full-reference guard passed: {full['guard_passed']}. Official prep and independent full-payload audit passed. It has not been uploaded, and no leaderboard improvement is claimed. Compact synthetic emission has substantial neutral-control DE bias. See [V10_STRATEGY.md](V10_STRATEGY.md) and [V10_SUBMISSION.md](V10_SUBMISSION.md). Reproduce with root script `run_v10_pipeline.py`.\n\n"
 if marker in text:text=text[:text.index(marker)]
 readme.write_text((text.rstrip()+'\n\n'+section).rstrip()+'\n',encoding='utf-8')
 if (A/'leaderboard_result.json').exists():
  from record_v10_result import document_publication
  document_publication()
 print('V10 documentation created from verified records',flush=True)
if __name__=='__main__':main()
