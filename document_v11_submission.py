"""Document V11 only after official preparation and independent full audit pass."""
import json,sys,platform
from importlib.metadata import version
import numpy as np,pandas as pd
from v11_methods import ROOT,A,OUT
import build_submission as build


def table(frame):
 cols=list(frame.columns);lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for _,row in frame.iterrows():lines.append('| '+' | '.join('undefined' if pd.isna(v) else f'{v:.6f}' if isinstance(v,float) else str(v) for v in row)+' |')
 return '\n'.join(lines)+'\n'


def main():
 def read(name):return json.loads((A/name).read_text(encoding='utf-8'))
 ready=read('submission_ready.json');audit=read('independent_validation.json');selection=read('count_selection.json');confirmation=read('confirmation_result.json');model=read('reference_model.json')
 tests=read('scientific_test_result.json');model_audit=read('model_validation.json')
 assert tests['passed'] and tests['tests']==4 and model_audit['passed']
 assert tests['test_script_sha256']==build.sha256(ROOT/'test_v11_submission.py') and tests['methods_sha256']==build.sha256(ROOT/'v11_methods.py')
 assert tests['output_sha256']==build.sha256(A/'scientific_test_output.txt')
 assert model_audit['model_sha256']==build.sha256(A/'reference_model.json')
 assert ready['official_prep_passed'] and audit['passed'] and audit['all_depths_preserve_sampled_controls']
 assert ready['sha256']==audit['package_sha256'] and ready['model_sha256']==audit['model_sha256']==build.sha256(A/'reference_model.json')
 assert ready['bytes']==(OUT/'prediction.vcc').stat().st_size and confirmation['confirmed_improvement']
 assert model['candidate']==selection['candidate']==confirmation['candidate']
 for stage in ['development','confirmation']:
  meta=read(f'{stage}_evaluation.json')
  for name,key in [(f'{stage}_scores.csv','scores_sha256'),(f'{stage}_panels.csv','panels_sha256'),(f'{stage}_summary.csv','summary_sha256')]:assert build.sha256(A/name)==meta[key]
  assert meta['inputs']['methods']==build.sha256(ROOT/'v11_methods.py') and meta['inputs']['script']==build.sha256(ROOT/'evaluate_v11_submission.py')
 assert selection['methods_sha256']==build.sha256(ROOT/'v11_methods.py') and selection['selection_script_sha256']==build.sha256(ROOT/'select_v11_submission.py')
 assert selection['plan_sha256']==build.sha256(A/'validation_plan.json') and confirmation['selection_sha256']==build.sha256(A/'count_selection.json')
 for name in ['development','confirmation']:
  power=read(f'{name}_power_validation.json')
  assert power['completed'] and power['methods_sha256']==build.sha256(ROOT/'v11_methods.py')
  assert power['script_sha256']==build.sha256(ROOT/'audit_v11_cell_power.py')
  assert power['scores_sha256']==build.sha256(A/f'{name}_power_scores.csv') and power['summary_sha256']==build.sha256(A/f'{name}_power_summary.csv')
 assert selection['development_power_validation_sha256']==build.sha256(A/'development_power_validation.json')
 assert confirmation['power_validation_sha256']==build.sha256(A/'confirmation_power_validation.json')
 neutral=read('decoder_neutral_manifest.json')
 assert neutral['scores_sha256']==build.sha256(A/'decoder_neutral_scores.csv') and neutral['summary_sha256']==build.sha256(A/'decoder_neutral_summary.csv')
 assert neutral['script_sha256']==build.sha256(ROOT/'benchmark_v11_neutral.py') and neutral['methods_sha256']==build.sha256(ROOT/'v11_methods.py')
 stage=read('memory_safe_stage.json');assert stage['source_sha256']==audit['prediction_sha256']
 assert ready['packaging_adapter_sha256']==stage['adapter_sha256']==build.sha256(ROOT/'package_v11_memory_safe.py')
 assert ready['memory_safe_stage_sha256']==build.sha256(A/'memory_safe_stage.json')
 fields=['reference_mode','candidate','pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
 text=f"# V11 verified submission\n\nSelected **{model['candidate']}**, setting `{json.dumps(model['setting'],sort_keys=True)}`. V10 remains the best published project submission at 0.116754121 until another official result beats it. V11 has not been uploaded.\n\n"
 text+='See [V11_STRATEGY.md](V11_STRATEGY.md) for the hypothesis, source review, code, frozen candidate grid and reproduction.\n\n'
 public=pd.read_csv(A/'confirmation_summary.csv').set_index(['reference_mode','candidate'])
 baseline=public.loc[('full_independent','v10')];candidate=public.loc[('full_independent',model['candidate'])]
 text+=f"Compared with V10 using full independent controls, public fidelity improves {baseline.fidelity:.6f} to {candidate.fidelity:.6f}, reach improves {baseline.reach:.6f} to {candidate.reach:.6f}, and DE fold-change error falls {baseline.nmae:.6f} to {candidate.nmae:.6f}. The common-control expression-error ratio worsens {100*(candidate.common_expression_ratio/baseline.common_expression_ratio-1):.2f}% ({baseline.common_expression_ratio:.6f} to {candidate.common_expression_ratio:.6f}; lower is better). These are public proxy measures, with the 400-cell SLIRP reach tradeoff and neutral-control bias detailed below.\n\n"
 for name in ['development','confirmation']:
  text+=f'## Public {name} results\n\n'+table(pd.read_csv(A/f'{name}_summary.csv')[fields])+'\n'
 text+=f"Frozen candidate confirmation guard passed: **{confirmation['guard_passed']}**. Weighted paired six-metric proxy change versus V10: **{confirmation['weighted_utility_difference']:+.6f}**. Public proxies are not hidden anchored scores. HCT/HEK count-confirmation targets were previously untested; H1 and public source summaries remain development evidence. No candidate was selected or retuned using confirmation.\n\n"
 rows=[]
 for mode,details in confirmation['details'].items():
  for family in ['H1','HCT116','HEK293T']:
   b=details['family_fidelity']['v10'][family];v=details['family_fidelity'][model['candidate']][family]
   rows.append({'reference_mode':mode,'family':family,'v10_fidelity':b,'v11_fidelity':v,'change':v-b})
 text+='## Fidelity by held-out family\n\n'+table(pd.DataFrame(rows))+'\n'
 text+='## 400-distinct-treated-cell H1 power check\n\n'+table(pd.read_csv(A/'confirmation_power_summary.csv'))+'\n'
 power_summary=pd.read_csv(A/'confirmation_power_summary.csv').set_index(['reference_mode','candidate'])
 power_base=power_summary.loc[('full_independent','v10')];power_new=power_summary.loc[('full_independent',model['candidate'])]
 text+=f'The power guard passed without duplicating cells. On SLIRP with full independent controls, fidelity changed {power_new.fidelity-power_base.fidelity:+.6f}, NMAE changed {power_new.nmae-power_base.nmae:+.6f}, but reach fell from {power_base.reach:.6f} to {power_new.reach:.6f}. The predefined power guard protects fidelity and NMAE, not reach. This target-specific tradeoff remains. Only H1 has enough observed treated cells here; this is reused public biology and does not remove the HCT/HEK or hidden-context power gap.\n\n'
 text+='## Neutral population bias\n\n'+table(pd.read_csv(A/'decoder_neutral_summary.csv'))+'\n'
 text+='Compact synthetic cells remain biased neutral populations. Reduced false calls do not mean faithful biological cell variance, and neither public guard passage nor this audit guarantees a leaderboard gain.\n\n'
 text+='## Production response changes\n\n'+table(pd.DataFrame(model['contexts']))+'\n'
 production=pd.read_csv(A/'prediction_group_statistics.csv')
 assert len(production)==900 and production.groupby(['context','target_gene']).size().eq(1).all()
 constraint_summary=production.groupby('context').agg(groups=('target_gene','size'),nonzero_counts=('nnz','sum'),max_mean_probability_L1=('mean_probability_L1','max'),max_bulk_probability_L1=('bulk_probability_L1','max')).reset_index()
 text+='## Production population constraints\n\n'+table(constraint_summary)+'\n'
 text+='L1 distances compare requested versus realized mean cell proportions and pooled proportions. All 900 groups meet the predefined limits; these constraints do not establish realistic cell variance.\n\n'
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 summary=pd.read_csv(A/'confirmation_summary.csv').set_index(['reference_mode','candidate']);fig,axes=plt.subplots(1,3,figsize=(12,3.6))
 for ax,metric,label in zip(axes,['fidelity','reach','nmae'],['Direction fidelity (higher better)','Direction reach (higher better)','DE fold-change error (lower better)']):
  for i,mode in enumerate(['sampled400','full_independent']):
   ax.bar(i-.16,summary.loc[(mode,'v10'),metric],width=.3,color='#687782',label='V10' if i==0 else None)
   ax.bar(i+.16,summary.loc[(mode,model['candidate']),metric],width=.3,color='#176b99',label='V11' if i==0 else None)
  ax.set_xticks([0,1],['400 controls','Full independent']);ax.set_title(label,fontsize=10);ax.spines[['top','right']].set_visible(False)
 axes[0].legend(frameon=False,fontsize=8);fig.suptitle('V11 count-confirmation: public raw proxies, not leaderboard scores',fontsize=11)
 fig.tight_layout();fig.savefig(A/'validation_comparison.png',dpi=180);plt.close(fig)
 text+=f"## Scientific and integrity checks\n\nFour focused scientific tests passed. The model audit verified {model_audit['evaluation_target_seed_groups']} target/seed splits and recomputed sampling for {model_audit['power400_groups']} power-check groups with 400 distinct treated cells. All production response arrays were independently recomputed exactly.\n\n"
 text+='![Public count-confirmation comparison](experiments/v11_consensus_coverage/analysis/validation_comparison.png)\n\n'
 text+=f"## Upload artifact\n\n- File: `experiments/v11_consensus_coverage/prediction.vcc`\n- Bytes: {ready['bytes']:,}\n- Shape: 360,000 cells x 18,533 genes; 900 groups x 400 cells.\n- Official prep, every packaged count/index/pointer/axis, all control provenance and depths: passed.\n- Exact replay: {audit['exact_replayed_groups']} groups, {audit['exact_replayed_gene_counts']:,} gene counts.\n- SHA-256: `{ready['sha256']}`\n- Uploaded: false. Leaderboard performance unknown.\n\n"
 slash=chr(92);cwd=slash.join(['L:','promec','Animesh','Download','Virtual']);exe=slash.join(['.','.competition_env','Scripts','python.exe']);file=slash.join(['.','experiments','v11_consensus_coverage','prediction.vcc'])
 text+=f'## Submit from Windows Command Prompt\n\n```bat\ncd /d {cwd}\n{exe} -m vcc submit {file} -m "Consensus-supported coverage and dispersion v11" --wait\n```\n\n'
 text+='If waiting is interrupted after upload, follow the existing entry with `vcc status <entry-id> --wait` instead of creating a duplicate.\n'
 (ROOT/'V11_SUBMISSION.md').write_text(text.rstrip()+'\n',encoding='utf-8')
 path=ROOT/'README.md';current=path.read_text(encoding='utf-8');marker='<!-- V11_VERIFIED_SUBMISSION -->'
 if marker in current:current=current[:current.index(marker)]
 current=current.rstrip()+f"\n\n{marker}\n\n## V11 upload candidate\n\nV11 is ready at `experiments/v11_consensus_coverage/prediction.vcc`, selected `{model['candidate']}`. Frozen count-confirmation passed; public weighted utility change **{confirmation['weighted_utility_difference']:+.6f}** versus V10. Official prep and independent full-payload/replay/depth audits passed. It has not been uploaded; V10 remains our best published entry. Synthetic neutral-control bias remains. See [V11_STRATEGY.md](V11_STRATEGY.md) and [V11_SUBMISSION.md](V11_SUBMISSION.md). Reproduce with root `run_v11_pipeline.py --resume`.\n"
 path.write_text(current,encoding='utf-8')
 packages=['numpy','scipy','pandas','h5py','anndata','scikit-learn','threadpoolctl','zstandard','vcc-cli','matplotlib']
 build.write_json(A/'runtime_versions.json',{'python':sys.version,'platform':platform.platform(),'packages':{p:version(p) for p in packages}})
 build.write_json(A/'documentation_validation.json',{'passed':True,'package_sha256':ready['sha256'],
  'script_sha256':build.sha256(ROOT/'document_v11_submission.py'),'submission_document_sha256':build.sha256(ROOT/'V11_SUBMISSION.md'),
  'selection_sha256':build.sha256(A/'count_selection.json'),'confirmation_sha256':build.sha256(A/'confirmation_result.json'),
  'all_corrected_evaluation_neutral_and_packaging_links_match':True,'scope':'Final linked provenance and report; no hidden performance claim.'})
 print('V11 verified documentation created',flush=True)
if __name__=='__main__':main()
